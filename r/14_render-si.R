source("r/header.R")

quarto::quarto_render("ms/si.qmd")

si = read_lines("ms/si.tex")
start_longtable = which(str_detect(si, "begin\\{longtable\\}"))
end_longtable = which(str_detect(si, "end\\{longtable\\}"))

# Change column widths
colwidths = c(
">{\\raggedright\\arraybackslash}p{(\\columnwidth - 8\\tabcolsep) * \\real{0.29}}",
">{\\raggedright\\arraybackslash}p{(\\columnwidth - 8\\tabcolsep) * \\real{0.08}}",
">{\\raggedright\\arraybackslash}p{(\\columnwidth - 8\\tabcolsep) * \\real{0.1}}",
">{\\raggedright\\arraybackslash}p{(\\columnwidth - 8\\tabcolsep) * \\real{0.20}}",
">{\\raggedright\\arraybackslash}p{(\\columnwidth - 8\\tabcolsep) * \\real{0.33}}@{}}"
)

si[start_longtable + seq(5)] <- colwidths

c(si[seq(start_longtable) - 1],
  "\\begin{landscape}",
  si[start_longtable:end_longtable],
  "\\end{landscape}",
  si[(end_longtable + 1):length(si)]) |>
  write_lines("ms/si.tex")

system('cd ms; xelatex si.tex; xelatex si.tex')
