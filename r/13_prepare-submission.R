# Move lines around so okina formats properly and then recompile
ms = readr::read_lines("ms/ms.tex")
x = which(ms %in% c(
  "\\usepackage{newunicodechar,graphicx}",
  "\\DeclareRobustCommand{\\okina}{\\raisebox{\\dimexpr\\fontcharht\\font`A-\\height}{\\scalebox{0.8}{`}}}",
  "\\newunicodechar{ʻ}{\\okina}"
))

readr::write_lines(c(ms[1:6], ms[x], ms[7:(min(x) - 1)], ms[(max(x) + 1): length(ms)]), "ms/ms.tex")
system("cd ms; xelatex ms.tex")

# Copy and rename files for submission
file.copy("ms/ms.pdf", "aobp1/ms.pdf", overwrite = TRUE)
file.copy("ms/si.pdf", "aobp1/si.pdf", overwrite = TRUE)
# file.copy("figures/study-system.pdf", "aobp1/fig1.pdf", overwrite = TRUE)
# file.copy("figures/habitat-aa.pdf", "aobp1/fig2.pdf", overwrite = TRUE)
# file.copy("figures/traits-aa.pdf", "aobp1/fig3.pdf", overwrite = TRUE)
# file.copy("figures/licor.pdf", "aobp1/figS1.pdf", overwrite = TRUE)
# file.copy("figures/pp-licor.pdf", "aobp1/figS2.pdf", overwrite = TRUE)
# file.copy("figures/habitat-Ags.pdf", "aobp1/figS3.pdf", overwrite = TRUE)
# file.copy("figures/habitat-gmaxratio.pdf", "aobp1/figS4.pdf", overwrite = TRUE)

# zip("aobp1/dryad.zip", list.files("dryad", full.names = TRUE))

system("cd aobp1; git clone git@github.com:cdmuir/stomata-spacing.git")

zip("aobp1/github.zip", list.files("aobp1/stomata-spacing", full.names = TRUE,
                                  recursive = TRUE))

system("rm -r aobp1/stomata-spacing")
