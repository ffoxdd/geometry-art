# Craftcloud

[Craftcloud](https://craftcloud3d.com/) is a 3D-printing marketplace: one
upload is quoted by many manufacturing partners, each price is binding, and
the order goes to whichever offer is picked. This note records what it
accepts, what drives its prices, and how to design a part to be cheap on it.
Facts are as read from its help centre, material guide and public API on
2026-09-25.

## What it accepts

- **One body per file.** A file with several separate bodies -- a
  "multi-shell" file -- is rejected. Each part is uploaded as its own file.
  Parts joined by tabs or sprues into one solid count as one body
  ([Multiple parts in one 3D file](https://support.craftcloud3d.com/en/articles/80-multiple-parts-in-one-3d-file)).
- **No customer supports.** Manufacturers orient the part and generate
  supports for their own machines, and charge the removal as labour
  ([Supports](https://support.craftcloud3d.com/en/articles/17-what-are-supports-what-do-i-need-to-know-about-them)).
- **Layer height** is 0.1 to 0.2 mm for plastics and down to 0.05 mm for
  high-detail resin. FDM infill is chosen at checkout; SLS and MJF parts
  print solid.
- **Raised or engraved text** needs at least 0.5 mm of height and width.
- **Watertight, manifold meshes.** Upload warnings flag open or
  self-intersecting geometry
  ([Upload warnings](https://support.craftcloud3d.com/en/articles/59-i-get-a-warning-when-uploading-my-files)).

The [Material Guide](https://craftcloud3d.com/en/material-guide) gives a
minimum wall and minimum detail per process
([Design Guide](https://support.craftcloud3d.com/en/articles/32-design-guide)):

| Process | Minimum wall | Minimum detail |
|---|---|---|
| FDM | 1.0 mm | 0.8 mm |
| SLS | 0.8 mm | 0.3 mm |
| MJF | 0.8 mm | 0.25 mm |
| SLA, MSLA, DLP, LCD | 0.6 mm | 0.25 mm |
| PolyJet | 1.0 mm | 0.25 mm |
| Binder jetting | 0.8 mm | 0.5 mm |

A feature at the minimum is printable but fragile. Sizing about a quarter
above it leaves margin for handling and shipping.

## What drives price

A quote depends on the part's **volume**, its **bounding box** and its
**surface area**, on the material and machine time of the process, and on
**labour** such as support removal and finishing
([Understanding Prices](https://support.craftcloud3d.com/en/articles/24-understanding-prices-on-craftcloud)).
Each vendor also has a **minimum production price**, returned alongside every
quote in the API.

What an order actually costs, from quotes to a US address:

- **An order's total is one vendor's charge plus its shipping.** The charge
  is the sum of the part prices, raised to the vendor's minimum production
  price when it falls short. Minimums run from nothing to nearly $200;
  the cheapest vendors sit around $25 to $40. Standard shipping runs about
  $9 to $19.
- **For small plastic parts the fixed charges are the price.** A flat
  15 cm FDM lattice of about 6 cm³ quotes near $3 as a part, and 32 of them
  near $16, so a small order pays mostly minimum and shipping.
- **So fill the order up to the minimum.** Past a handful of parts, more
  parts cost almost nothing until their sum passes the minimum.
- **Separate parts in one order cost no more than a joined sheet.** Four
  separate swatches quote below one sheet of four, because the vendor minimum
  applies to the order, not to each part.
- **The cheapest materials** are the plain FDM thermoplastics: PLA, PLA+,
  PETG and PCTG. SLS nylon is the cheapest powder process, at several times
  the price for the same part.
- **The website's material list totals every uploaded model**, in that
  material, across the whole upload. Per-model prices show only after
  selecting a material, or in the API.

How the terms bear on a design:

- **Volume** is material. Thin members cost in proportion to cross-section,
  so a lattice's volume scales with bar width times thickness.
- **Bounding box** is machine space. Powder-bed processes (SLS, MJF) pack
  many customers' parts into one build, so a part with a large box and
  little material still occupies the space it encloses.
- **Surface area** tracks finishing and, for resin, peel and exposure work.
- **Per-order fixed charges** favour putting everything wanted into one
  order.

## Layout principles

- **Lay thin parts flat and side by side.** A flat sheet wastes only the gaps
  between parts. Stacking pays an air gap above every layer, and the gap must
  be wide enough for the shop to clear powder or resin, so for parts about
  1 mm thick a stack roughly doubles the box.
- **Keep parts separate** unless a design needs them joined. The minimum is
  charged per order, so joining saves nothing and adds warping risk and bed
  constraints.
- **Watch the bed size.** Common FDM beds are about 250 mm, large ones about
  350 mm. An HP MJF 5200 builds in 380 × 284 × 380 mm. A part bigger than a
  bed can only be quoted by the shops with bigger machines.
- **Large thin flat parts warp.** They are the worst case for FDM and SLS
  shrinkage, which argues for several medium sheets over one maximal one.

## API

The quoting API is described at
<https://api.craftcloud3d.com/api-docs.json> (v5).

- `POST /model` uploads a mesh and returns a model id.
- `POST /v5/price` takes a currency, a country code and a list of models
  with quantity and scale, optionally narrowed to material configurations or
  vendors, and returns a price id.
- `GET /v5/price/{priceId}` polls until `allComplete`, returning quotes per
  vendor and material with shipping options and each vendor's minimum
  production price.

The price endpoints answer a plain HTTP client, but the model upload sits
behind a Cloudflare bot challenge. A model uploaded through the website keeps
its id in the browser's local storage, under the session's backend models,
and that id can be priced from anywhere. The `all3dp` Ruby gem
(<https://github.com/cults/all3dp>) wraps a different endpoint: it hands
Craftcloud public URLs of mesh files and gets back a configuration link for
a person to open. It returns no prices.
