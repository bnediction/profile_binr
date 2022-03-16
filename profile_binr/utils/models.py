"""
    Models for bifurcation trajectory reconstruction.
"""

from .customobjs import ObjDict

CELL_FATE = ObjDict(
    {
        "network": ObjDict(
            {  # modèle réduit
                "TNF": "TNF",
                "FAS": "FAS",
                "RIP1": "!C8 & (TNF | FAS)",
                "NFkB": "(cIAP & RIP1) & !C3",
                "C8": "(TNF | FAS | C3) & !NFkB",
                "cIAP": "(NFkB | cIAP) & !MOMP",
                "ATP": "!MPT",
                "C3": "ATP & MOMP & !NFkB",
                "ROS": "!NFkB & (RIP1 | MPT)",
                "MOMP": "MPT | (C8 & !NFkB)",
                "MPT": "ROS & !NFkB",
                "Apoptosis": "C3",
                "NonACD": "!ATP",
                "Survival": "NFkB",
            }
        ),
        "init_active": ["cIAP", "ATP", "TNF"],  # nœuds actifs initialement
    }
)

EARLY_HEMATOPOIESIS = ObjDict(
    {
        "network": ObjDict(
            {
                "Egr1": "Gata2 & Junb",
                "Junb": "Egr1 | Myc",
                "Bclaf1": "Myc",
                "Myc": "Cebpa & Bclaf1",
                "Fli1": "Junb | (Gata1 & !Klf1)",
                "Gata2": "(Gata2 & !Gata1 & !Zfpm1) | (Egr1 & !Gata1 & !Zfpm1 & !Spi1)",
                "Spi1": "(Spi1 & !Gata1) | (Cebpa & !Gata1 & !Gata2)",
                "Cebpa": "(Gata2 & !Ikzf1) | (Spi1 & !Ikzf1)",
                "Gata1": "Fli1 | (Gata2 & !Spi1) | (Gata1 & !Ikzf1 & !Spi1)",
                "Klf1": "Gata1 & !Fli1",
                "Tal1": "Gata1 & !Spi1",
                "Ikzf1": "Gata2",
                "Zfpm1": "Gata1",
                "CDK46CycD": "Bclaf1 | Myc",
                "CIPKIP": "Junb",
            }
        ),
        "init_active": [],  # TODO : complete this
    }
)

CORE_REGULATION_1 = ObjDict(
    {
        "network": ObjDict(
            {
                "G1": "TF1",
                "G2": "TF1",
                "G3": "TF1",
                "G4": "TF2",
                "G5": "TF2",
                "G6": "TF2",
                "TF1": "TF3 & !TF2",
                "TF2": "TF3 & !TF1",
                "TF3": "TF4",
                "TF4": "TF5",
                "TF5": "TF6",
                "TF6": 1,
            }
        ),
        "init_active": ["TF6"],
    }
)
