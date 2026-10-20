# Token and call usage audit

Project folder: `/workspace/nucleus/work/projects/p-3b5b7d6e8e314ceaa070504180c1227b`
Repository: `repo` (branch `main`, 61 commits)
Session id seen in job records: `can-you-make-a-wfg8`
Audit performed: read-only. No discovery run, swarm or long computation was started.
No energy or CO2 conversion is made anywhere in this document, and none should be inferred from it.

## Date range covered

| Boundary | Timestamp (UTC) | What it is |
|---|---|---|
| Earliest timestamp visible | 2026-10-19 17:20 (minute precision) | project directory mtime |
| Earliest *event* timestamp | 2026-10-19T17:26:36Z | first retained shell job start |
| Earliest commit | 2026-10-19T19:16:03Z | `43f763c`, "Core entropy of the real quadratic family" |
| Latest job event | 2026-10-20T17:50:42Z | last retained shell job end |
| Latest timestamp visible | 2026-10-20T23:34:07Z | `cdc37a7`, last commit |

Span: 2026-10-19T17:20Z to 2026-10-20T23:34:07Z = **30 h 14 min** (1.26 days).

Normalising to a rate: every quantity the service measured for this project is zero
(see table), so the rate of measured spend is 0 tokens/day, 0 model calls/day,
0 USD/day over this span. That is a statement about *what is recorded*, not a claim
that no model work happened — see "What is not measurable here".

## Sources

| Source | Date range | Model calls | Tokens in | Tokens out | USD | Confidence | Evidence |
|---|---|---|---|---|---|---|---|
| Discovery runs, `discovery_list` (whole project) | none — no run exists | 0 | 0 tokens | 0 tokens | 0 USD | measured-by-service | `discovery_list` returned `{"runs": []}`; no `next` cursor, so the listing is exhausted in one page |
| Discovery per-run detail, `discovery_status` / `discovery_problems` / `discovery_tree` | n/a | not available | not available | not available | not available | not available | No run id exists to query; a probe of `discovery_problems` was refused with `no such run` |
| Swarms this session, `ultramath_status` / `ultramath_records` | none | 0 | not available | not available | 0 USD | measured-by-service | Both tools refused: `this session has started no swarm` |
| Swarms from other sessions, on disk | none found | not available | not available | not available | not available | counted-from-files | No `.ultramath/`, `swarms/` or `runs/` directory anywhere under the project folder (`ls`, `find`) |
| Run journals / transcripts / manifests | none found | not available | not available | not available | not available | not available | `rg -l -i 'prompt_tokens\|model_calls\|completion_tokens\|total_tokens\|tokens_spent\|\busd\b'` over the whole project folder (hidden files included, `.git` excluded): no matches |
| Shell job records, `.nucleus/jobs/` | 2026-10-19T17:26:36Z – 2026-10-20T17:50:42Z | not available | not available | not available | not available | counted-from-files (bytes/seconds only) | 8 log files, **2,967 bytes** total, **238.6 s** summed wall clock; fields present are path, size, state, exit code, tail — no token or call fields |
| Conversation turns and their tool round trips | n/a | not available | not available | not available | not available | not available | No transcript, message store or turn counter is readable from this project |

Column meaning: "Tokens in" = prompt/input tokens, "Tokens out" = completion/output tokens.
`0` above is a measured emptiness (the source exists and reports nothing for this project),
not an estimate; `not available` means the source carries no such field.

### Shell job records, itemised (bytes only)

| Job id | Started (UTC) | Ended (UTC) | State | Exit | Log bytes |
|---|---|---|---|---|---|
| `jmvfj0mf2` | 2026-10-19T17:26:36 | 2026-10-19T17:26:42 | exited | 1 | 210 |
| `jmvfj6xkl` | 2026-10-19T17:31:31 | 2026-10-19T17:31:47 | exited | 0 | 1,180 |
| `jmvgz5h1n` | 2026-10-20T17:46:03 | 2026-10-20T17:46:03 | exited | 1 | 161 |
| `jmvgz5r4b` | 2026-10-20T17:46:16 | 2026-10-20T17:47:37 | exited | 0 | 195 |
| `jmvgz8wmi` | 2026-10-20T17:48:43 | 2026-10-20T17:49:11 | stopped by agent | — | 0 |
| `jmvgz9moo` | 2026-10-20T17:49:17 | 2026-10-20T17:50:42 | exited | 0 | 1,196 |
| `jmvgz9xoh` | 2026-10-20T17:49:31 | 2026-10-20T17:50:37 | exited | 0 | 12 |
| `jmvgzajlm` | 2026-10-20T17:49:59 | 2026-10-20T17:50:15 | exited | 0 | 13 |

Total 2,967 bytes. These logs are the only machine-side record of shell output that
survives; they bound *retained* tool-output volume at the low end only, because most
tool output in a session is returned inline and never written to a job log.

Other byte counts observed, for scale (not usage):

- uploaded artefact `repo/attachments/can-you-make-a-wfg8/88MT:IteratedInterval.pdf`: 5,378,491 bytes
- repo working tree including `.git`: 47 MiB; tracked files: 16

## Tool calls and turns evidenced

- **Tool calls with a surviving record: 8.** Each retained shell job corresponds to one
  `shell_exec` call with `background: true` (7 exited, 1 stopped by the agent). That is a
  hard lower bound on tool calls in this project and nothing more.
- **Assistant turns: not derivable.** No source counts turns. 61 commits witness 61
  assistant actions that produced a commit, but a commit is not a turn and not a tool call,
  and the number of turns behind each commit is unknown.
- **Model (LLM) calls: not derivable** from any source. Nothing in the project records them.
- This audit itself used **13 read-only tool calls** to gather the material above (listings,
  `rg`, `du`, `git log`, timestamp arithmetic), plus 2 more to write and verify this file —
  15 in total. It started no run and consumed no prover budget.

## Model identity

**Not determinable from anything in this project.** There are no run records, no worker
manifests, no config file and no journal that names a serving or prover model — the
search above for run/usage fields returned nothing, and neither `discovery_list` (empty)
nor the `ultramath_*` tools (no swarm) expose a model string for this project. The
experiment-policy manifests under `/opt/nucleus/share/policies` were deliberately **not**
consulted: they name the model that built the exploration policy, which is a different
artefact from the serving model.

## What is not measurable here

1. **Ordinary conversational turns and their tool-call round trips appear in no source
   above.** The audit sees exactly three kinds of artefact: service-side run/ledger data
   (empty for this project), background-job logs (8 files, 2,967 bytes, no usage fields),
   and files the work produced (`repo`). The prompt, the completion, the tool call
   arguments and the tool results of every interactive turn are absent from all three.
   This is the dominant gap, and it is not a small remainder.
2. **Size of that blind spot — what can and cannot be bounded.** It cannot be bounded in
   tokens, and no honest number can be produced for it: there is no per-turn ledger, no
   message store and no character count of the transcript. What can be said concretely is
   that one input artefact is known in exact size (a 5,378,491-byte PDF was uploaded into
   the session), that 61 commits and 8 background jobs were produced inside a 30-hour span,
   and that a session of that shape necessarily carries a great many turns whose input
   alone includes the repository sources and that uploaded document on every turn it was
   in context. Any specific token figure for the conversational channel would be an
   invention; this audit does not supply one.
3. **Discovery and swarm spend for other sessions of this project.** The tools are
   per-session or service-scoped, and empty. Filesystem evidence shows no `.ultramath/`,
   `swarms/` or `runs/` directory, so on-disk records for other sessions do not exist under
   this project folder. If such runs exist anywhere else, they are not visible from here
   and are reported as `not available`, not as zero.
4. **Model identity** — see the previous section.
5. **Tool-output volume beyond job logs.** Inline tool results (file reads, `rg` output,
   Lean or Python output returned in-call) are not retained anywhere, so total tool-output
   bytes are unbounded above; 2,967 bytes is only the retained background-job share.
6. **Wall-clock time is not call volume.** 238.6 s of summed background-job time and the
   30 h 14 min project span describe elapsed time, and neither maps to model calls without
   a per-call ledger that does not exist.

## Summary

Every source that actually measures spend reports nothing for this project: `discovery_list`
is empty, no swarm exists in this session or on disk, and no journal, transcript or manifest
carries a token or call field. What exists is developmental output — 61 commits in 30 h 14 min
produced with 8 retained background jobs totalling 2,967 bytes of log — plus a 5,378,491-byte
uploaded PDF. The measured-but-empty sources are trustworthy: they are the service's own
ledgers, and they were paged to exhaustion. The conversational channel, which is almost
certainly where nearly all model calls and tokens in this project went, is recorded in no
source reachable from here, and it is the single unknown that dominates the total.

---

### Method / reproducibility

Commands run, in order (all read-only except the write of this file):

- `discovery_list`; `ultramath_status`; `ultramath_records`; `shell_jobs list`
- `discovery_problems` with an invented run name, to confirm the refusal mode
- `ls -la`, `find .nucleus -type f -printf`, `du -sb/-bc`, `du -sh` on the project folder
  and `.nucleus/{clones,uploads,bundles,trash}` (all empty)
- `rg --files --hidden` and `rg -l -i '<usage field pattern>'` over the whole project folder,
  `.git` excluded
- `git -C repo log --format='%H %ad %s' --date=iso-strict`, `git rev-list --count HEAD`
- `python3` (standard library only) to convert the epoch-millisecond job timestamps to UTC
  and to sum job wall clock
- write of this file with `fs_create`

Nothing was downloaded, no network call was made, and no run was started.

### Scope note

`usage-audit.md` sits in the project folder, which is **not** a git repository — `repo/` is
the only repository here. This file therefore does not appear in the `git diff` review view,
by the deliberate choice not to initialise a repository the person did not ask for.
