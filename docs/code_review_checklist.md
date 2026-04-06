# Code Review Checklist — Scientific Computing Team

> **Team:** 2–5 engineers | **Stack:** Python, Fortran, C, C++, Shell, TensorFlow, HAPI, ExoMol
> **Context:** Most code is AI-generated. Reviews must validate correctness at the code level — not just at the plan/report level — and produce an audit trail that can answer a PI's detailed questions.

---

## 0. AI-Generated Code Validation

- [ ] **Provenance is documented.** Commit message or PR description states which parts were AI-generated vs. hand-written.
- [ ] **No hallucinated APIs or functions.** Every library call, function name, and method signature actually exists in the version we use. *(Spot-check at least 3 calls against docs.)*
- [ ] **No plausible-but-wrong logic.** AI code often "looks right" but implements a subtly different algorithm. Verify core logic against the specification or paper, not just by reading it.
- [ ] **No dead or vestigial code.** Unused imports, commented-out blocks, or redundant branches are removed.
- [ ] **No hardcoded assumptions.** Check for magic numbers, hardcoded paths, or assumptions about array shapes/sizes that should be parameterized.
- [ ] **Tests were written or updated.** AI-generated code is not considered reviewed until its behavior is verified by a test.
- [ ] **Reviewer actually ran the code.** At least one reviewer executed the code (or its tests) locally.

---

## 1. Correctness

- [ ] **Algorithm matches specification.** Implementation matches the referenced paper, equation, or design doc.
- [ ] **Edge cases are handled.** Empty inputs, zero-length arrays, NaN/Inf values, single-element cases, boundary conditions.
- [ ] **Error handling is explicit.** Functions fail loudly (raise/assert) rather than silently returning wrong results.
- [ ] **Units and dimensions are consistent.** Physical quantities carry clear units. Conversions are explicit.
- [ ] **Floating-point precision.** Appropriate types used. No naive summation of large float arrays where it matters.
- [ ] **Off-by-one errors.** Array indexing is correct across language boundaries.
- [ ] **Concurrency and race conditions.** Shared state is protected, file writes don't collide.

---

## 2. Security

- [ ] **No secrets in code.** API keys, tokens, passwords not committed.
- [ ] **No unsafe shell calls.** `subprocess.run()` with argument lists, not `shell=True`.
- [ ] **File paths are validated.** No path traversal vulnerabilities.
- [ ] **Dependencies are pinned and audited.** Exact versions in requirements.
- [ ] **Input data is validated.** External files checked for format, size, content.
- [ ] **Permissions are appropriate.** Output files not world-writable. Temp files cleaned up.

---

## 3. Performance

- [ ] **Computational complexity is acceptable.** No O(n^2) or worse on large spectral datasets without justification.
- [ ] **Memory usage is bounded.** Large arrays not duplicated unnecessarily.
- [ ] **I/O is not a bottleneck.** Binary formats preferred over repeated text parsing for large data.
- [ ] **Vectorization is used where possible.** NumPy operations replace explicit Python loops.
- [ ] **No redundant computation.** Expensive calculations cached or precomputed.
- [ ] **Parallelism is effective.** Speedup measured and documented.

---

## 4. Data Pipeline Integrity

- [ ] **Input file format is validated on read.** Column counts, header keywords, expected ranges checked.
- [ ] **Output file format matches downstream expectations.**
- [ ] **Intermediate files are documented.** What they contain, when deletable, what depends on them.
- [ ] **File I/O is atomic or recoverable.** Long-running writes use temp files + rename.
- [ ] **Data types are preserved across boundaries.** No silent float32→float64, no integer truncation.

---

## 5. Reproducibility

- [ ] **Random seeds are set and logged.**
- [ ] **Environment is captured.** Package versions, compiler versions, OS info recorded.
- [ ] **Results can be regenerated.** Same inputs → identical output.
- [ ] **Version of code is tagged.** Commit hash recorded alongside output.
- [ ] **External data sources are versioned.** ExoMol/HITRAN versions noted.

---

## 6. Maintainability

- [ ] **Naming is clear and domain-appropriate.** `temperature_K`, `cross_section_cm2`, not `t`, `xs`.
- [ ] **Functions are focused.** Each function does one thing.
- [ ] **Comments explain *why*, not *what*.**
- [ ] **Fortran/C interop is documented.** Argument types, array ordering, memory ownership.
- [ ] **No copy-paste duplication.** Repeated blocks extracted into functions.
- [ ] **The code is navigable.** File and module structure is logical.

---

## 7. Process & Audit Trail

- [ ] **The PR links to its plan/task.**
- [ ] **The commit history tells a story.** Logical units with explanatory messages.
- [ ] **Review comments are preserved.** In the PR, not in ephemeral chat.
- [ ] **Deviations from the plan are documented.**
- [ ] **Known limitations are stated.**
- [ ] **Results are linked.** Report references exact code version and input data.

---

*Last updated: 2026-04-02 | Review and revise quarterly or after any incident.*
