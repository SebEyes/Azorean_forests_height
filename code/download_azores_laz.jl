#!/usr/bin/env julia
#
# download_azores_laz.jl
#
# Download all Azores LiDAR point-cloud tiles (DGT collection ACORES-LAZ)
# from https://cdd.dgterritorio.gov.pt/dgt-fe/mapa?collection=ACORES-LAZ
#
# The portal has no documented API, but its frontend uses a STAC API and a
# Keycloak-protected download endpoint. This script:
#   1. enumerates every tile with the STAC search API (pagination included),
#   2. logs in through the same Keycloak OIDC flow used by the website
#      (or reuses a browser `connect.sid` cookie),
#   3. downloads each .laz file and verifies it against the LAS header / STAC
#      point count (`file:size` is only advisory — DGT metadata can be stale).
#
# Requires: HTTP, JSON3 (see code/00-julia_env_initialisation.jl).
#
# Examples
# --------
#   # list what would be downloaded (no login needed)
#   julia code/download_azores_laz.jl --list-only
#
#   # download everything (login with your DGT account)
#   julia -t 4 code/download_azores_laz.jl \
#       --username USER --password PASS --outdir data/azores_laz
#
#   # instead of credentials, reuse a cookie copied from your browser
#   julia code/download_azores_laz.jl --cookie '<connect.sid value>' \
#       --outdir data/azores_laz
#
#   # one island only, first 5 tiles, preview
#   julia code/download_azores_laz.jl --island SMG --limit 5 --dry-run

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))
Pkg.instantiate()

using HTTP
using JSON3
using Dates
using Printf: @sprintf
using Base.Threads: @threads, threadid, nthreads

const SITE        = "https://cdd.dgterritorio.gov.pt"
const API         = "$SITE/dgt-be/v1"
const SEARCH_URL  = "$API/search"
const DOWNLOAD_URL = "$API/download"
const AUTH_START  = "$SITE/auth/login"
const CLIENT_ID   = "aai-oidc-dgt"   # not sent explicitly: the server starts the PKCE flow
const COLLECTION  = "ACORES-LAZ"
const PAGE_SIZE   = 1000             # STAC limit; API caps at 1000
const OUTFILE_EXT = ".laz"

# Network tuning -------------------------------------------------------------
# HTTP.jl's default `connect_timeout` is only 30 s and it also bounds the TLS
# handshake. The DGT download CDN can take considerably longer than that,
# which surfaces as HTTP.TimeoutError("tls_handshake", 30_000_000_000, ...).
# Use a generous per-attempt connect timeout (overridable with
# --connect-timeout), plus explicit limits for the other request phases.
const CONNECT_TIMEOUT         = Ref(120)    # s: TCP connect + TLS handshake (HTTP.jl requires an Int)
const RESPONSE_HEADER_TIMEOUT = 120.0       # s: waiting for response headers
const READ_IDLE_TIMEOUT       = 300.0       # s: max inactivity while streaming
const DOWNLOAD_ATTEMPTS       = 5           # attempts per tile (fresh link each)
const RETRY_BACKOFF           = 2.0         # s, doubled after each failed attempt
const LAS_MAGIC               = "LASF"      # LAS/LAZ files all start with this

# Island codes as they appear in item IDs (e.g. 2024_SMG_... = São Miguel).
const ISLAND_NAMES = Dict(
    "SMG" => "São Miguel",
    "SMA" => "Santa Maria",
    "TER" => "Terceira",
    "GRA" => "Graciosa",
    "SJO" => "São Jorge",
    "PIC" => "Pico",
    "FAI" => "Faial",
    "FLO" => "Flores",
    "COR" => "Corvo",
)

const ISLAND_CODES = Dict{String,String}()
for (code, name) in ISLAND_NAMES
    ISLAND_CODES[lowercase(name)] = code
    ISLAND_CODES[lowercase(replace(name, "ã" => "a", "é" => "e", "í" => "i",
                                        "ó" => "o", "ú" => "u", "ç" => "c"))] = code
    ISLAND_CODES[lowercase(code)] = code
end
ISLAND_CODES["sjg"] = "SJO"   # alias used in code/00-julia_env_initialisation.jl
for k in collect(keys(ISLAND_CODES))
    occursin(' ', k) && (ISLAND_CODES[replace(k, " " => "-")] = ISLAND_CODES[k])
end

# Cookies are collected while following redirects and then reused for API calls.
const COOKIES     = Dict{String,String}()
const COOKIE_LOCK = ReentrantLock()   # COOKIES is touched by every worker thread
const PRINT_LOCK  = ReentrantLock()
const STATS       = Dict{Symbol,Int}(:ok => 0, :skip => 0, :fail => 0)

# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

jget(o, key::AbstractString, default=nothing) =
    haskey(o, Symbol(key)) ? getproperty(o, Symbol(key)) : default

function hget(r, name::AbstractString)
    for (k, v) in r.headers
        lowercase(k) == lowercase(name) && return v
    end
    return nothing
end

html_unescape(s::AbstractString) =
    replace(String(s), "&amp;" => "&", "&quot;" => "\"", "&#39;" => "'",
                      "&lt;" => "<", "&gt;" => ">")

function store_cookies!(r)
    lock(COOKIE_LOCK) do
        for (k, v) in r.headers
            lowercase(k) == "set-cookie" || continue
            first = split(v, ';')[1]
            occursin('=', first) || continue
            name, value = split(first, '=', limit=2)
            COOKIES[String(name)] = String(value)
        end
    end
    return nothing
end

function cookie_header(url::AbstractString)
    host = HTTP.URIs.parse_uri(url).host
    occursin("dgterritorio.gov.pt", host) || return ""
    return lock(COOKIE_LOCK) do
        isempty(COOKIES) ? "" : join(["$(k)=$(v)" for (k, v) in COOKIES], "; ")
    end
end

function request(method::AbstractString, url::AbstractString;
                 headers=Pair{String,String}[], body=nothing,
                 response_stream=nothing, redirect=false, status_exception=false,
                 retry::Bool=true, retries::Int=4, closeimmediately::Bool=false)
    hs = Pair{String,String}[headers...]
    c = cookie_header(url)
    isempty(c) || push!(hs, "Cookie" => c)

    # `connect_timeout` covers both the TCP connect and the TLS handshake in
    # HTTP.jl v2. Without it the 30 s default is used for each attempt.
    kw = Dict{Symbol,Any}(:redirect => redirect,
                          :status_exception => status_exception,
                          :connect_timeout => CONNECT_TIMEOUT[],
                          :response_header_timeout => RESPONSE_HEADER_TIMEOUT,
                          :read_idle_timeout => READ_IDLE_TIMEOUT,
                          :retry => retry,
                          :retries => retries,
                          :closeimmediately => closeimmediately)
    body !== nothing && (kw[:body] = body)
    response_stream !== nothing && (kw[:response_stream] = response_stream)

    r = HTTP.request(method, url, hs; kw...)
    store_cookies!(r)
    return r
end

resolve_location(loc::AbstractString, base::AbstractString) =
    startswith(loc, "http") ? loc :
    startswith(loc, "/")    ? SITE * loc :
    replace(base, r"/[^/]*$" => "/") * loc

formencode(fields) =
    join(["$(HTTP.URIs.escapeuri(string(k)))=$(HTTP.URIs.escapeuri(string(v)))"
          for (k, v) in fields], "&")

# Follow redirects, returning the final page text, response and URL.
function get_page(url::AbstractString; maxhops::Int=10)
    cur = url
    for _ in 1:maxhops
        buf = IOBuffer()
        r = request("GET", cur; response_stream=buf)
        loc = hget(r, "Location")
        if 300 <= r.status < 400 && loc !== nothing
            cur = resolve_location(loc, cur)
            continue
        end
        return String(take!(buf)), r, cur
    end
    error("Too many redirects while fetching $url")
end

# ---------------------------------------------------------------------------
# Keycloak login (same browser flow as the portal; PKCE is handled by the
# portal itself: we only need to keep the cookies across redirects).
# ---------------------------------------------------------------------------

function parse_form_element(form::AbstractString)
    am = match(r"action=\"([^\"]+)\"", form)
    am === nothing && return nothing
    fields = Dict{String,String}()
    for im in eachmatch(r"<input[^>]*>", form)
        tag = im.match
        nm = match(r"name=\"([^\"]+)\"", tag)
        nm === nothing && continue
        tm = match(r"type=\"([^\"]+)\"", tag)
        typ = tm === nothing ? "" : lowercase(tm.captures[1])
        vm = match(r"value=\"([^\"]*)\"", tag)
        val = vm === nothing ? "" : html_unescape(vm.captures[1])
        name = html_unescape(nm.captures[1])
        if typ in ("hidden", "submit") || name in ("username", "password", "otp", "credentialId")
            fields[name] = val
        end
    end
    return html_unescape(am.captures[1]), fields
end

function parse_form(html::AbstractString, id::AbstractString)
    m = match(Regex("<form[^>]*id=\"$id\".*?</form>", "s"), html)
    m === nothing ? nothing : parse_form_element(m.match)
end

function parse_otp_form(html::AbstractString)
    for m in eachmatch(r"<form[^>]*>.*?</form>"s, html)
        occursin("name=\"otp\"", m.match) || continue
        parsed = parse_form_element(m.match)
        parsed !== nothing && return parsed
    end
    return nothing
end

function submit_form(action::AbstractString, fields::Dict{String,String})
    buf = IOBuffer()
    r = request("POST", action;
                headers=["Content-Type" => "application/x-www-form-urlencoded"],
                body=formencode(fields), response_stream=buf)

    cur = action
    for _ in 1:10
        loc = hget(r, "Location")
        if 300 <= r.status < 400 && loc !== nothing
            cur = resolve_location(loc, cur)
            buf = IOBuffer()
            r = request("GET", cur; response_stream=buf)
            continue
        end
        break
    end
    return String(take!(buf)), r, cur
end

function verify_login()
    r = request("GET", "$API/auth/user"; response_stream=IOBuffer())
    loc = hget(r, "Location")
    (r.status == 200 && (loc === nothing || !occursin("/auth/login", loc))) && return true

    # Fallback: the download route itself is behind the auth middleware.
    r2 = request("GET", "$DOWNLOAD_URL/__session_probe__"; response_stream=IOBuffer())
    loc2 = hget(r2, "Location")
    return loc2 === nothing || !occursin("/auth/login", loc2)
end

function login!(username::AbstractString, password::AbstractString; otp=nothing)
    html, _, _ = get_page(AUTH_START)

    parsed = parse_form(html, "kc-form-login")
    parsed === nothing && error("Could not find the Keycloak login form; the portal may have changed.")

    action, fields = parsed
    fields["username"] = username
    fields["password"] = password
    fields["login"] = get(fields, "login", "Entrar")
    haskey(fields, "credentialId") || (fields["credentialId"] = "")

    html, r, cur = submit_form(action, fields)

    # Two-step verification, if enabled on the account.
    if occursin("name=\"otp\"", html)
        parsed_otp = parse_otp_form(html)
        parsed_otp === nothing && error("OTP requested but no OTP form found.")
        otp === nothing && error("This account requires a 2FA code. Pass --otp 123456.")
        oact, ofields = parsed_otp
        ofields["otp"] = String(otp)
        html, r, cur = submit_form(oact, ofields)
    end

    if occursin("id=\"kc-form-login\"", html)
        error("Login failed: check --username/--password (Keycloak rejected the credentials).")
    end

    verify_login() || error("Login did not produce an authenticated portal session.")
    println("  ✓ Logged in as $username")
    return nothing
end

# ---------------------------------------------------------------------------
# STAC enumeration
# ---------------------------------------------------------------------------

island_code(id::AbstractString) =
    (m = match(r"^\d{4}_([A-Z]+)_", id)) === nothing ? "UNK" : m.captures[1]

function collect_items(; island=nothing, limit=nothing)
    body = Dict{String,Any}(
        "filter-lang" => "cql2-json",
        "filter" => Dict{String,Any}(
            "op" => "eq",
            "args" => Any[Dict("property" => "collection"), COLLECTION],
        ),
        "limit" => PAGE_SIZE,
    )

    items = Any[]
    while true
        r = request("POST", SEARCH_URL;
                    headers=["Content-Type" => "application/json"],
                    body=JSON3.write(body))
        r.status == 200 || error("STAC search returned HTTP $(r.status): $(String(r.body)[1:min(end,300)])")
        data = JSON3.read(String(r.body))
        append!(items, collect(data.features))

        token = nothing
        for l in data.links
            jget(l, "rel", "") == "next" || continue
            b = jget(l, "body", nothing)
            b === nothing && continue
            t = jget(b, "token", nothing)
            t !== nothing && (token = String(t))
        end

        token === nothing && break
        body["token"] = token
    end

    if island !== nothing
        code = get(ISLAND_CODES, lowercase(String(island)), uppercase(String(island)))
        filter!(i -> island_code(String(i.id)) == code, items)
    end
    limit !== nothing && resize!(items, min(length(items), limit))
    sort!(items; by=i -> String(i.id))
    return items
end

function print_summary(items)
    bycount = Dict{String,Int}()
    bysize  = Dict{String,Float64}()
    total   = 0.0
    for i in items
        code = island_code(String(i.id))
        sz = Int(jget(jget(i, "properties"), "file:size", 0))
        bycount[code] = get(bycount, code, 0) + 1
        bysize[code]  = get(bysize, code, 0.0) + sz
        total += sz
    end
    println("\nFound $(length(items)) tile(s), $(round(total/1e9, digits=1)) GB total:")
    for code in sort(collect(keys(bycount)))
        name = get(ISLAND_NAMES, code, code)
        println(@sprintf("  %3s  %-14s %5d files  %7.1f GB",
                         code, name, bycount[code], bysize[code]/1e9))
    end
    println()
end

# ---------------------------------------------------------------------------
# Downloading
# ---------------------------------------------------------------------------

# The portal's download links are single-use tokens bound to the session that
# performed the STAC search. A failed/partially consumed token therefore has to
# be replaced by a brand-new one before retrying.
function fresh_asset_href(item)
    id = String(item.id)
    refreshed = collect_items(; island=island_code(id))
    for it in refreshed
        String(it.id) == id || continue
        assets = jget(it, "assets")
        asset = jget(assets, "data", jget(assets, "Data", nothing))
        asset === nothing && continue
        return String(asset.href)
    end
    error("Could not refresh the download link for $id")
end

function asset_href(item)
    assets = jget(item, "assets")
    asset = jget(assets, "data", jget(assets, "Data", nothing))
    asset === nothing && error("No download asset found for $(String(item.id))")
    return String(asset.href)
end

"""
    las_header_info(path)

Read the authoritative metadata from a LAS/LAZ file header (no decompression).
The STAC `file:size` field can be stale by a few hundred bytes when the DGT
re-compresses a tile, so the embedded point count is used for validation instead.
"""
function las_header_info(path::AbstractString)
    open(path, "r") do io
        b = read(io, 512)
        length(b) >= 256 || error("Truncated LAS header in $(basename(path))")
        String(b[1:4]) == LAS_MAGIC ||
            error("Not a LAS/LAZ file (starts with $(repr(String(b[1:4]))))")

        major, minor = Int(b[25]), Int(b[26])
        (major == 1 && 0 <= minor <= 4) || error("Unsupported LAS version $major.$minor")

        # LAS header byte offsets (0-based here; Julia indices are +1).
        hdrsize = Int(ltoh(reinterpret(UInt16, b[95:96])[1]))    # 94
        offset  = Int(ltoh(reinterpret(UInt32, b[97:100])[1]))   # 96
        fmt     = b[105]                                          # 104
        legacy  = Int(ltoh(reinterpret(UInt32, b[108:111])[1]))  # 107
        # LAS 1.4 has separate 64-bit extended counts at byte 247.
        total = (major == 1 && minor >= 4) ?
            Int(ltoh(reinterpret(UInt64, b[248:255])[1])) : legacy

        return (version=(major, minor), format=Int(fmt & 0x3f),
                compressed=(fmt & 0x80) != 0, hdrsize=hdrsize,
                offset=offset, count=total)
    end
end

"""
    download_to_file(href, path)

Stream a GET response to `path`, following redirects manually.

HTTP.jl assigns `response_stream` to every request in its redirect chain, so
using `redirect=true` writes the intermediate `302 Found. Redirecting to …`
body into the output file before the actual LAZ payload. Resolving each hop
here keeps only the final 200 response on disk.
"""
function download_to_file(href::AbstractString, path::AbstractString)
    url = href
    for hop in 1:10
        rm(path; force=true)
        next = nothing
        open(path, "w") do io
            # A fresh connection per hop: HTTP.jl pools sockets, and a socket
            # reused after a streamed 302 response can already be closed.
            r = request("GET", url; response_stream=io, redirect=false, retry=false,
                        closeimmediately=true)
            if r.status in (301, 302, 303, 307, 308)
                loc = HTTP.header(r, "Location")
                isempty(loc) && error("Redirect without Location for $(basename(path))")
                next = resolve_location(loc, url)
            else
                r.status in (200, 206) ||
                    error("HTTP $(r.status) while downloading $(basename(path))")
            end
        end
        next === nothing && return nothing
        url = next
    end
    error("Too many redirects while downloading $(basename(path))")
end

function verify_download(tmp::AbstractString, id::AbstractString,
                         expected_size::Int, expected_count::Int)
    sz = filesize(tmp)
    sz == 0 && error("Empty download for $id")

    info = las_header_info(tmp)
    info.offset > 0 && info.offset < sz ||
        error("Truncated LAS/LAZ file for $id (size $sz, point offset $(info.offset))")
    info.count > 0 || error("LAS/LAZ file for $id reports zero points")
    expected_count > 0 && info.count != expected_count &&
        error("Point count mismatch for $id: file has $(info.count), STAC says $expected_count")

    # file:size is advisory: it occasionally differs from the served file by a
    # few hundred bytes. The point count above is the real integrity check.
    if expected_size > 0 && abs(sz - expected_size) > max(4096, ceil(Int, 0.001 * expected_size))
        @warn "Download size for $id differs from STAC metadata: got $sz, " *
              "expected $expected_size (point count $(info.count) validated)"
    end
    return nothing
end

function valid_existing(path::AbstractString, expected_size::Int, expected_count::Int)
    try
        info = las_header_info(path)
        filesize(path) > info.offset || return false
        (expected_count == 0 || info.count == expected_count) || return false
        # A truncated file with a valid header is still rejected using a very
        # conservative lower bound (>= 0.5 byte per compressed point).
        filesize(path) >= info.offset + max(1, info.count ÷ 2) || return false
        return true
    catch
        return false
    end
end

function download_item(item; outdir::AbstractString, force::Bool=false, dryrun::Bool=false)
    id   = String(item.id)
    props = jget(item, "properties")
    expected = Int(jget(props, "file:size", 0))
    expected_count = Int(jget(props, "pc:count", 0))
    path = joinpath(outdir, id * OUTFILE_EXT)

    if !force && isfile(path) && valid_existing(path, expected, expected_count)
        return :skip, expected
    end
    dryrun && return :would, expected

    href = asset_href(item)
    tmp = path * ".part.$(threadid())"
    last_err = nothing

    for attempt in 1:DOWNLOAD_ATTEMPTS
        # Fresh temporary file each time: a half-written file must never be
        # mistaken for a complete download, and must never be appended to.
        isfile(tmp) && rm(tmp; force=true)
        try
            if attempt > 1
                # Expired/consumed tokens need a new link; back off first so a
                # struggling CDN is not hammered.
                sleep(RETRY_BACKOFF * 2^(attempt - 2))
                href = fresh_asset_href(item)
            end

            # Do not probe with HEAD: it consumes the single-use token.
            # Redirects are followed manually so the 302 body is not written
            # into the file. Retries are handled by this loop so every attempt
            # gets a fresh link and a fresh temporary file.
            download_to_file(href, tmp)

            verify_download(tmp, id, expected, expected_count)
            mv(tmp, path; force=true)
            return :ok, expected
        catch e
            last_err = e
            isfile(tmp) && rm(tmp; force=true)
            if attempt < DOWNLOAD_ATTEMPTS
                @warn "Attempt $attempt/$DOWNLOAD_ATTEMPTS failed for $id: $e — retrying"
            end
        end
    end

    error("Failed to download $id after $DOWNLOAD_ATTEMPTS attempts: $last_err")
end

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

function parse_cli()
    opts = Dict{String,String}()
    flags = Set{String}()
    i = 1
    while i <= length(ARGS)
        a = ARGS[i]
        if startswith(a, "--")
            if occursin('=', a)
                k, v = split(a[3:end], '=', limit=2)
                opts[k] = v
            else
                k = a[3:end]
                if k in ("help", "list-only", "dry-run", "force")
                    push!(flags, k)
                else
                    i += 1
                    i <= length(ARGS) || error("Missing value for --$k")
                    opts[k] = ARGS[i]
                end
            end
        end
        i += 1
    end
    return opts, flags
end

function help()
    println("""
    Usage: julia [-t N] code/download_azores_laz.jl [options]

    Options:
      --username USER     DGT portal username/email (or env DGT_USERNAME)
      --password PASS     DGT portal password (or env DGT_PASSWORD)
      --otp CODE          Two-factor code (if your account uses 2FA)
      --cookie VALUE      Existing `connect.sid` cookie copied from the browser
      --outdir DIR        Destination directory (default: data/azores_laz)
      --island CODE       Limit to one island: SMG SMA TER GRA SJO PIC FAI FLO COR
                          (names such as "sao-miguel" are also accepted)
      --limit N           Download at most N tiles
      --workers N         Parallel downloads (also set julia -t N), default 1
      --connect-timeout N Connect/TLS-handshake timeout in seconds (default 120);
                          raise this if downloads fail with tls_handshake timeouts
      --list-only         Only list tiles/sizes; no login or downloads
      --dry-run           Show which files would be downloaded
      --force             Re-download files even if they already exist
      --help              Show this help
    """)
end

function main()
    opts, flags = parse_cli()
    "help" in flags && (help(); return)

    outdir   = get(opts, "outdir", joinpath("data", "azores_laz"))
    island   = get(opts, "island", nothing)
    limit    = haskey(opts, "limit") ? parse(Int, opts["limit"]) : nothing
    workers  = haskey(opts, "workers") ? parse(Int, opts["workers"]) : 1
    haskey(opts, "connect-timeout") && (CONNECT_TIMEOUT[] = parse(Int, opts["connect-timeout"]))
    CONNECT_TIMEOUT[] > 0 || error("--connect-timeout must be a positive number of seconds")
    dryrun   = "dry-run" in flags
    listonly = "list-only" in flags
    force    = "force" in flags

    println("Querying STAC catalogue for $COLLECTION ...")
    items = collect_items(; island=island, limit=limit)
    isempty(items) && (println("No tiles matched."); return)
    print_summary(items)

    listonly && return

    if dryrun
        println("Files that would be written to $(abspath(outdir)):")
        for item in items
            println("  ", joinpath(outdir, String(item.id) * OUTFILE_EXT))
        end
        return
    end

    username = get(opts, "username", get(ENV, "DGT_USERNAME", ""))
    password = get(opts, "password", get(ENV, "DGT_PASSWORD", ""))

    if haskey(opts, "cookie")
        COOKIES["connect.sid"] = opts["cookie"]
        verify_login() || error("The supplied connect.sid is not authenticated (or has expired).")
        println("  ✓ Using browser session cookie")
    elseif !isempty(username) && !isempty(password)
        login!(username, password; otp=get(opts, "otp", get(ENV, "DGT_OTP", nothing)))
    else
        error("Provide --username/--password (or DGT_USERNAME/DGT_PASSWORD) or --cookie (or use --list-only).")
    end

    # Download hashes are bound to the session that performed the STAC search.
    # Login regenerates the session id, so the hashes collected above are stale.
    println("  Refreshing download links for the authenticated session ...")
    items = collect_items(; island=island, limit=limit)

    mkpath(outdir)
    println("Downloading to $(abspath(outdir)) with $workers worker(s)...\n")

    failed = String[]
    if workers <= 1
        for (n, item) in enumerate(items)
            id = String(item.id)
            try
                status, bytes = download_item(item; outdir=outdir, force=force, dryrun=dryrun)
                lock(PRINT_LOCK) do
                    get!(STATS, status, 0); STATS[status] += 1
                    status === :skip || println("[$n/$(length(items))] $status  $id  " *
                                                @sprintf("%.0f MB", bytes/1e6))
                end
            catch e
                push!(failed, id)
                lock(PRINT_LOCK) do
                    STATS[:fail] += 1
                    println(stderr, "[$n/$(length(items))] FAIL  $id : $e")
                end
            end
        end
    else
        worker_count = min(workers, nthreads())
        worker_count < workers && @warn "julia was started with $(nthreads()) thread(s); using $worker_count."
        @threads for n in eachindex(items)
            item = items[n]
            id = String(item.id)
            try
                status, bytes = download_item(item; outdir=outdir, force=force, dryrun=dryrun)
                lock(PRINT_LOCK) do
                    get!(STATS, status, 0); STATS[status] += 1
                    status === :skip || println("[$n/$(length(items))] $status  $id  " *
                                                @sprintf("%.0f MB", bytes/1e6))
                end
            catch e
                lock(PRINT_LOCK) do
                    STATS[:fail] += 1
                    push!(failed, id)
                    println(stderr, "[$n/$(length(items))] FAIL  $id : $e")
                end
            end
        end
    end

    println("\nDone: $(get(STATS,:ok,0)) downloaded, $(get(STATS,:skip,0)) already present, " *
            "$(get(STATS,:fail,0)) failed.")
    isempty(failed) || println("Failed IDs (retry the script, completed files are skipped):\n  " *
                               join(failed, "\n  "))
end

# Only run the CLI when executed directly; when `include`d, expose the helper
# functions (STAC enumeration, login, download_item, ...) for other scripts.
(abspath(PROGRAM_FILE) == @__FILE__) && main()

