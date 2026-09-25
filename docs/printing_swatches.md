# Printing swatches

A swatch is a square plane tessellation, thickened by `skeletonize` into a
lattice of bars with a full-width frame along its walls. Swatches exist to be
experimented on -- covered, coated, cast over -- so the print is a substrate:
the cheapest material that holds its shape is the right one, and cost is
what the layout optimises. They are ordered through Craftcloud, whose rules
and pricing are in [craftcloud.md](craftcloud.md).

The working target is a swatch about 15 cm on a side with cells whose clear
opening is at least about 1 cm. On a constant density that is roughly 180
cells; a gradient or radial field needs more points to keep its sparsest
cells at that opening.

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

## Quote experiment

Each row is one upload, quoted in the cheapest material that meets its
process minimum, in one currency and destination:

| Model | What it isolates | Cheapest quote |
|---|---|---|
| swatch, 1.2 mm bars | the baseline | |
| swatch, 1.2 mm bars, quantity 2 and 4 | per-part against per-order charges | |
| swatch, 1.0 and 1.6 mm bars | volume and surface area at a near-equal box | |
| solid plate, same box | volume at an equal box | |
| sheet 2 × 1 | two swatches as one part | |
| sheet 3 × 1 | a part longer than common FDM beds | |
| sheet 2 × 2 | a part wider than common FDM beds | |

Comparing a sheet of n swatches with n copies of one swatch prices grouping
directly; comparing the bar sizes and the solid plate separates the volume
term from the box term.
