# CG__250425_10x — Raw Data

Raw 10x Genomics FASTQs for the **L364** experimens (4 captures, A–D) plus two top-up runs for L364C.
Each capture has three library types: gene expression (GEX), BCR V(D)J, and cell surface protein (CSP: ADT/HTO).

---

## Layout

```
Raw_Data/
├── Run_<SeqRun>_<Capture>[_TopUpSequencing_NN]/
│   ├── GEX/        # Gene expression          (*_GEX_*)
│   ├── BCR/        # BCR V(D)J                (*_VDJ_*)
│   └── ADT_HTO/    # Antibody / hashtag CSP   (*_CSP_*)
```

## Runs

| Folder | Seq run | Capture | Library date | Lanes | GEX | BCR | ADT_HTO |
|---|---|---|---|---|---|---|---|
| `Run_S1081_L364A` | S1081 | L364A | 21-02-2025 | L001 | ✓ | ✓ | ✓ |
| `Run_S1146_L364B` | S1146 | L364B | 08-07-2025 | L001–L002 | ✓ | ✓ | ✓ |
| `Run_S1146_L364C` | S1146 | L364C | 15-07-2025 | L001–L002 | ✓ | ✓ | ✓ |
| `Run_S1146_L364D` | S1146 | L364D | 25-07-2025 | L001–L002 | ✓ | ✓ | ✓ |
| `Run_S1152_L364C_TopUpSequencing_01` | S1152 | L364C | 15-07-2025 | L001 | – | – | ✓ |
| `Run_S1161_L364C_TopUpSequencing_02` | S1161 | L364C | 15-07-2025 | L001–L002 | – | – | ✓ |

## File naming

```
[<SeqRun>_]<Capture>_<DDMMYYYY>_<Library>_S<n>_L00<lane>_<Read>_001.fastq.gz
```

| Field | Values |
|---|---|
| `SeqRun` | prefix only on top-up runs (`S1152`, `S1161`) |
| `Capture` | `L364A`–`L364D` |
| `Library` | `GEX`, `VDJ`, `CSP` |
| `Read` | `I1`, `I2` (indices), `R1` (barcode+UMI), `R2` (insert) |

## Notes for processing

- **L364C CSP** is spread over three runs (S1146, S1152, S1161). Pass all three `ADT_HTO/` folders as comma-separated `fastqs` for the same feature library in `cellranger multi`; sample prefixes differ (`L364C_15072025_CSP` vs `S1152_L364C_15072025_CSP` / `S1161_L364C_15072025_CSP`), so list each prefix explicitly.