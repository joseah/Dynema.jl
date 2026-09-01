# ---------------------------------------------------------------------------- #
#                                cli_output.jl                                 #
# ---------------------------------------------------------------------------- #
#
# Small formatting helpers for consistent, readable progress output across
# the CLI scripts: labeled section headers with indented bullet points
# instead of a flat stream of unlabeled println's, plus a `with_tee_log`
# helper that duplicates everything printed during a run to a log file
# (alongside the exact command invoked) while still showing it on the
# console as normal. `include()`d before vcf_geno.jl/mtx_expr.jl by every
# entry-point script that uses them (dynema_map.jl, dynema_extract_geno.jl),
# since those shared files call `section`/`bullet` too.

using Dates

"""
    section(msg)

Prints a top-level progress header, e.g.:

    ==> Loading metadata
"""
section(msg::AbstractString) = println("\n==> ", msg)

"""
    bullet(msg; indent=1)

Prints an indented bullet point under the current `section`, e.g.:

    - file: metadata.tsv

Nest further (e.g. for a breakdown under one bullet) with `indent=2`, etc.
"""
bullet(msg::AbstractString; indent::Int = 1) = println("    "^indent, "- ", msg)

"""
    elapsed_str(t0) -> String

Formats the elapsed wall-clock time since `t0 = time()` as e.g. `"2.35s"`,
for a trailing `bullet("done in \$(elapsed_str(t0))")`.
"""
elapsed_str(t0::Real) = string(round(time() - t0, digits = 2), "s")

"""
    shell_quote(s) -> String

Single-quotes `s` (escaping any embedded single quotes) if it contains
whitespace or other shell-special characters, otherwise returns it
unchanged. Used to log the invoking command in a form that can be
copy-pasted and re-run as-is -- important here since paths under this
project routinely contain spaces (e.g. ".../Partners HealthCare Dropbox/...").
"""
function shell_quote(s::AbstractString)
    occursin(r"[^A-Za-z0-9_./=:,+-]", s) || return s
    return "'" * replace(s, "'" => raw"'\''") * "'"
end

"""
    strip_ansi(s) -> String

Removes ANSI/VT100 escape sequences (SGR color/style codes like the ones
Dynema's own `show` method prints, e.g. `Base.printstyled`-generated
`\\e[93;1m`) from `s`. Terminals interpret these as colors, but a plain text
log file just shows the raw escape bytes -- so `with_tee_log` uses this to
keep the log human-readable while still letting the live console output be
colored.
"""
strip_ansi(s::AbstractString) = replace(s, r"\e\[[0-9;]*[A-Za-z]" => "")

"""
    with_tee_log(f, logpath; command)

Runs `f()` with stdout and stderr duplicated to both the terminal (exactly
as they'd print without this wrapper, ANSI colors included) and a fresh log
file at `logpath` (with ANSI codes stripped via `strip_ansi`, since they're
just clutter in a plain text file), so a complete transcript of one run --
progress sections/bullets, `@warn` messages, and the final results table --
is saved to disk. `command` is written as a header line first, so the log
is self-contained: it records exactly how the run was invoked, not just
what it printed.

Implementation note: `redirect_stdout()`/`redirect_stderr()` (no arguments)
create OS-level pipes and make Julia's `stdout`/`stderr` (and therefore
`println`, `@warn`, `show`, etc.) write into them; a background task reads
each pipe line-by-line and re-emits it to both the original terminal stream
and the log file. The original streams are restored in a `finally` block so
this cleans up properly even if `f()` errors (the error itself, printed by
Julia after `with_tee_log` returns/rethrows, will land on the terminal but
not in the log -- see the log's own last lines / the terminal for that).
"""
function with_tee_log(f, logpath::AbstractString; command::AbstractString)

    log_io = open(logpath, "w")
    println(log_io, "# Command:  ", command)
    println(log_io, "# Started:  ", Dates.now())
    flush(log_io)

    orig_stdout, orig_stderr = stdout, stderr
    rd_out, wr_out = redirect_stdout()
    rd_err, wr_err = redirect_stderr()

    tee(rd, orig) = @async try
        for line in eachline(rd)
            println(orig, line)
            println(log_io, strip_ansi(line))
            flush(log_io)
        end
    catch
        # Reader task dies naturally once the write end is closed below;
        # nothing to surface if that happens mid-iteration.
    end
    t_out = tee(rd_out, orig_stdout)
    t_err = tee(rd_err, orig_stderr)

    try
        f()
    finally
        redirect_stdout(orig_stdout)
        redirect_stderr(orig_stderr)
        close(wr_out)
        close(wr_err)
        wait(t_out)
        wait(t_err)
        println(log_io, "# Finished: ", Dates.now())
        close(log_io)
    end

end
