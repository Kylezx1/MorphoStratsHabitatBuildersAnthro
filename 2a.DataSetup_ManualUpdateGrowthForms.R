#Manually set growth form for some LIT species


LITDat$value[which(LITDat$species %in%
                     c(
                       "Acropora arborescent",
                       "Acropora formosa",
                       "Acropora robustea"
                     ))] <- "Arborescent"

LITDat$value[which(LITDat$species %in% c(
  "Acropora tabular"
))] <- "Tabular"

LITDat$value[which(LITDat$species %in% c(
  "Astreopora myriopthalma",
  "Faviidae",
  "Favites sp.",
  "Favites sp",
  "Favia pallida",
  "Favities magnistellata",
  "Favia stelligera",
  "Favites abdita?",
  "Lobophyllia hemphrichii",
  "Porites sp.",
  "Porites massive",
  "Porites",
  "Porites sp",
  "Platygyra sp.",
  "Goniastrea sp.",
  "Goniastrea sp",
  "Montastraea curta",
  "Dipsastraea truncatus",
  "Dipsastraea mathaii",
  "Dipsastraea matthai",
  "Dipsastraea rotundata",
  "Symphyllia radians",
  "Coscinarea columna"
))] <- "Massive"

LITDat$value[which(LITDat$species %in% c(
  "Acropora fat_digitifera",
  "Acropora digitifera_fat",
  "Acropora sp. Fat dig",
  "Acropora tall_monticulosa"
))] <- "Digitate"

LITDat$value[which(LITDat$species %in% c(
  "Acropora spathulata"
))] <- "Corymbose"

LITDat$value[which(LITDat$species %in% c(
  "Acropora palifera",
  "Isopora",
  "Isoporan",
  "Isopora pallida"

))] <- "Columnar"

LITDat$value[which(LITDat$species %in% c(
  "Echinopora gemmacae"
))] <- "Laminar"

LITDat$value[which(LITDat$species %in% c(
  "Stylophora pistulata"
))] <- "BranchingClosed"


LITDat$value[which(LITDat$species %in% c(
  "Leptastrae transversa",
  "Platygyra yaeyamaensis"
))] <- "Submassive"

