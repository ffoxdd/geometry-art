# Printing swatches

A swatch is a square cut from a larger tessellation, thickened by
`skeletonize` into a lattice of bars with a frame around it. Swatches exist to be
experimented on -- covered, coated, cast over -- so the print is a substrate:
the cheapest material that holds its shape is the right one, and cost is
what the layout optimises. They are ordered through Craftcloud, whose rules
and pricing are in [craftcloud.md](craftcloud.md).

The working target is a swatch about 15 cm on a side with cells whose clear
opening is at least about 1 cm. On a constant density that is roughly 180
cells; a gradient or radial field needs more points to keep its sparsest
cells at that opening.

## Windows

A tessellation's cells deform against a wall, so a swatch is not the whole
domain but a window onto its middle: `skeletonize --window <w>` keeps a
centred square `w` wide, lays it flat, and runs a frame of
`--frame-width` (the bar width unless set) around it, outside. The cells'
openings run up to the window's edge, so the window is exactly the region
that shows, and the part is the window plus two frame widths across.

A torus has no walls at all, which makes it the best source: its window
shows the tessellation as if it went on forever. The domain only has to
hold the window and frame; a margin of a cell or two is plenty. A 6 cm
window, cut from a torus 7 cm on a side, sized in millimetres:

```
tessellate -g torus --width 70 --height 70 -f constant -p 55 --snapshot swatch
skeletonize swatch.json -w 1.2 -t 1.2 --window 60
```

## Bar size

Bars about a quarter above the process minimum: 1.2 mm for FDM, 1.0 mm for
SLS and MJF, 0.8 mm for resin.

A 180-cell constant-density swatch, 150 mm on a side, square bars:

| Bars | Volume | Surface area | Bounding box | Box volume |
|---|---|---|---|---|
| 1.0 mm | 3.98 cm³ | 156 cm² | 151.0 × 151.0 × 1.0 mm | 22.8 cm³ |
| 1.2 mm | 5.68 cm³ | 185 cm² | 151.2 × 151.2 × 1.2 mm | 27.4 cm³ |
| 1.6 mm | 9.93 cm³ | 240 cm² | 151.6 × 151.6 × 1.6 mm | 36.8 cm³ |

Volume scales with bar width times thickness and with total edge length;
edge length grows with the square root of the cell count at a fixed size.
The lattice fills about a fifth of its bounding box.

## Sheets

Swatches lie in a grid, frame to frame, with two tabs across each 3 mm gap
that are cut after printing. Tabs add under a percent of material. Joined
sheets, 1.2 mm bars:

| Sheet | Volume | Bounding box | Box volume |
|---|---|---|---|
| 2 × 1 | 11.38 cm³ | 305.4 × 151.2 × 1.2 mm | 55.4 cm³ |
| 3 × 1 | 17.08 cm³ | 459.6 × 151.2 × 1.2 mm | 83.4 cm³ |
| 2 × 2 | 22.80 cm³ | 305.4 × 305.4 × 1.2 mm | 111.9 cm³ |

## What swatches cost

Cheapest order total per count, to a US address, 1.2 mm bars in the cheapest
material on offer, quoted 2026-09-25. Every size keeps its cell openings
near 1 cm: 28 cells at 6 cm, 80 at 10 cm, 180 at 15 cm.

| Swatches | 6 cm | 10 cm | 15 cm |
|---|---|---|---|
| 1 | $18.43 | $18.43 | $36.43 |
| 2 | $22.16 | $22.16 | $37.83 |
| 4 | $29.61 | $29.61 | $42.25 |
| 8 | $38.70 | $42.94 | $44.50 |
| 16 | $43.99 | $52.29 | $59.03 |
| 32 | $54.32 | $59.03 | $65.08 |
| 64 | $59.03 | $71.40 | $78.49 |

Repeated requests for the same order differ by up to about a fifth, as
vendors join or drop out of a quote, so these are indicative.

Those counts are copies of one design, which carry a quantity discount.
Distinct designs do not: each separate model is handled and priced on its
own, so distinct swatches belong on a shared sheet.

Sixteen distinct 6 cm swatches, 1.2 mm bars:

| Ordered as | Best total |
|---|---|
| sixteen separate models | $55.67 |
| one tabbed 4 × 4 sheet, 25 cm square | $20.48 |

Larger sheets of 6 cm swatches, each one model:

| Sheet | Side | Vendors quoting | Best total | Per swatch |
|---|---|---|---|---|
| 4 × 4 | 25 cm | 36 | $20.48 | $1.28 |
| 5 × 5 | 32 cm | 30 | $36.52 | $1.46 |
| 6 × 6 | 38 cm | 13 | $48.65 | $1.35 |
| 8 × 8 | 51 cm | 8 | $64.45 | $1.01 |

A 4 × 4 sheet is the working unit: it fits common FDM beds, most vendors
quote it, and it is among the cheapest per swatch. More swatches go on more
4 × 4 sheets in the same order; two copies of one sheet quoted at $26.26.

Part prices alone, cheapest vendor:

| Part | Price |
|---|---|
| 15 cm swatch, 1.0 mm bars | $2.56 |
| 15 cm swatch, 1.2 mm bars | $2.91 |
| 15 cm swatch, 1.6 mm bars | $3.82 |
| solid plate, same box | $4.51 |
| sheet 2 × 1 | $4.14 |
| sheet 3 × 1 | $5.35 |
| sheet 2 × 2 | $5.11 |
