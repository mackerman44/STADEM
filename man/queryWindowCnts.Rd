% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/queryWindowCnts.R
\name{queryWindowCnts}
\alias{queryWindowCnts}
\title{DART Window Counts}
\source{
\url{https://www.cbr.washington.edu/dart}
}
\usage{
queryWindowCnts(
  dam = c("LWG", "WFF", "BON", "TDA", "JDA", "MCN", "IHR", "LMN", "LGS", "PRO", "ROZ",
    "PRD", "WAN", "RIS", "TUM", "RRH", "WEL", "ZOS"),
  spp_code = c("fc", "fk", "fb", "fs", "fsw", "fa", "fcj", "fkj", "fbj", "fsj", "fl",
    "ft"),
  start_date = NULL,
  end_date = NULL
)
}
\arguments{
\item{dam}{the dam code for the dam you wish to query for window counts. Possible codes are: WFF (Willamette Falls), BON (Bonneville), TDA (The Dalles), JDA (John Day), MCN (McNary), IHR (Ice Harbor), LMN (Lower Monumental), LGS (Little Goose), LWG (Lower Granite), PRO (Prosser), ROZ (Roza), PRD (Priest Rapids), WAN (Wanapum), RIS (Rock Island), TUM (Tumwater), RRH (Rocky Reach), WEL (Wells), ZOS (Zosel)}

\item{spp_code}{species code(s) to query window counts for. Possible codes are: fc (Chinook), fk (Coho), fb (Sockeye), fs (Steelhead), fsw (Wild Steelhead), fa (Shad), fcj (Jack Chinook), fkj (Jack Coho), fbj (Jack Sockeye), fsj (Jack Steelhead), fl (Lamprey), ft (Bull Trout}

\item{start_date}{character vector of date (\code{YYYYMMDD}) when query should start}

\item{end_date}{character vector of date (\code{YYYYMMDD}) when query should end}
}
\description{
Query and download window counts at a dam by species and year using DART
}
\examples{
queryWindowCnts(spawn_yr = 2015)
}
\author{
Kevin See
}
