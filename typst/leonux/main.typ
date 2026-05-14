#import "@preview/leonux:1.1.0": *

// Setup the format and enter personal information
#show: setup.with(
	ratio: "16-9",
	primary: rgb("137C24"), // This is the default color, if primary: none
	title: "Study of stability post Thermal Quench States",
	subtitle: "JOREK simulations",
	date: "December 15, 2025",
	author: "Román García Guill",
	institute: "IPP",
)

// slide number only counts slide and content,references (both based on slide)

// Create titlepage
#titlepage()

// Create Table of contents (shows sections, not slides)
#content(title: "Contents")

// Create a section
#section(title: "Study Parameters")

// Create a slide
#slide(title: "Physics Parameters Selection")[
	- ZKPerp
	- Some citation @ref of Leonux #footnote("Leonux - minimalistic typst slides")
	// Create a block (e.g. for a definition)
	#my-block(title: "Definition: Law of large numbers")[
		$ forall epsilon > 0: lim_(n -> infinity) P(|r_n - p| <= epsilon) = 1 $
		Or in words: For large $n$ applies: $p approx r_n$.
	]

	This text will be shown on the first and second subslide. \
	#show: later // Make following content appear on the next subslide
	This text will only be shown on the second subslide.
]

// Create bibliography

#references(title:"References")[
	#bibliography("bibliography.bib")
]
