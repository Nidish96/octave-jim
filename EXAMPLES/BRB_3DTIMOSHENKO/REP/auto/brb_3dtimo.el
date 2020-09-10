(TeX-add-style-hook
 "brb_3dtimo"
 (lambda ()
   (TeX-run-style-hooks
    "latex2e"
    "beamertmd"
    "beamertmd10"
    "graphicx"
    "subcaption")
   (LaTeX-add-labels
    "fig:xfresp"
    "fig:yfresp"
    "fig:zfresp"))
 :latex)

