templates = {
    "Electrophilic addition into ene": {
        "scope": [
            {"name": "chlorine", "[R1X]": "[Cl&H1&+0:3]", "[R2X]": "[Cl&H0&+0:3]"},
            {"name": "bromine", "[R1X]": "[Br&H1&+0:3]", "[R2X]": "[Br&H0&+0:3]"},
            {"name": "iodine", "[R1X]": "[I&H1&+0:3]", "[R2X]": "[I&H0&+0:3]"},
            {"name": "fluorine", "[R1X]": "[F&H1&+0:3]", "[R2X]": "[F&H0&+0:3]"},
        ],
        "mechanism": [
            {
                "template": "[C;H2,H1,H0;+0:1]=[C&+0:2].[R1X]>>[C;H3,H2,H1+0:1]-[C&+0:2]-[R2X]",
                "description": "Electrophilic addition",
            }
        ],
    },
    # "Michael Additions (intramolecular from -ium)": {
    #     "scope": [
    #         {
    #             "name": "amine/phosphonic acid into iminium/sulfonium/oxonium",
    #             "[R1X]": "[N,P&H1;+0:5]",
    #             "[R2X]": "[N,P&H1;+1:5]",
    #             "[R3X]": "[N,S,O;+1:4]",
    #             "[R4X]": "[N,S,O;+0:4]",
    #         },
    #         {
    #             "name": "carbanion/hydroxide/alkoxide/thiol into iminium/sulfonium/oxonium",
    #             "[R1X]": "[C,O,S;!$(O[N+]);-1:5]",
    #             "[R2X]": "[C,O,S;!$(O[N+]);+0:5]",
    #             "[R3X]": "[N,S,O;+1:4]",
    #             "[R4X]": "[N,S,O;+0:4]",
    #         },
    #     ],
    #     "mechanism": [
    #         {
    #             "template": "[C&+0:1]=[R3X][*:8][*:7][R1X]>>[C&+0:1]1-[R4X][*:8][*:7][R2X]1",
    #             "description": "5-mem nuc addition into -ium",
    #         },
    #         {
    #             "template": "[C&+0:1]=[R3X][*:8][*:7][*:6][R1X]>>[C&+0:1]1-[R4X][*:8][*:7][*:6][R2X]1",
    #             "description": "6-mem nuc addition into -ium",
    #         },
    #         {
    #             "template": "[C&+0:1]=[R3X][*:8][*:7][*:6][*:9][R1X]>>[C&+0:1]1-[R4X][*:8][*:7][*:6][*:9][R2X]1",
    #             "description": "7-mem nuc addition into -ium",
    #         },
    #         {
    #             "template": "[C&+0:1]=[R3X][*:8][*:7][*:6][*:9][*:10][R1X]>>[C&+0:1]1-[R4X][*:8][*:7][*:6][*:9][*:10][R2X]1",
    #             "description": "8-mem nuc addition into -ium",
    #         },
    #     ],
    # },
    "Michael Additions": {
        "scope": [
            {
                "name": "amine/phosphonic acid into iminium/sulfonium/oxonium",
                "[R1X]": "[N&!$(NC=O),P&H1;+0:5]",
                "[R2X]": "[N&!$(NC=O),P&H1;+1:5]",
                "[R3X]": "[N,S,O;+1:4]",
                "[R4X]": "[N,S,O;+0:4]",
            },
            {
                "name": "carbanion/hydroxide/alkoxide/thiol into iminium/sulfonium/oxonium",
                "[R1X]": "[C,O,S;!$(O[N+]);-1:5]",
                "[R2X]": "[C,O,S;!$(O[N+]);+0:5]",
                "[R3X]": "[N,S,O;+1:4]",
                "[R4X]": "[N,S,O;+0:4]",
            },
            {
                "name": "amine/phosphonic acid into carbonyl/thial",
                "[R1X]": "[N&!$(NC=O),P&H1;+0:5]",
                "[R2X]": "[N&!$(NC=O),P&H1;+1:5]",
                "[R3X]": "[O,S;+0:4]",
                "[R4X]": "[O,S;-1:4]",
            },
            {
                "name": "carbanion/hydroxide/alkoxide/thiol into carbonyl/thial",
                "[R1X]": "[C,O,S;!$(O[N+]);-1:5]",
                "[R2X]": "[C,O,S;!$(O[N+]);+0:5]",
                "[R3X]": "[O,S;+0:4]",
                "[R4X]": "[O,S;-1:4]",
            },
        ],
        "mechanism": [
            {
                "template": "[R1X].[C;!$(C[*;-1]);!$(C(=O)O);+0:1]=[R3X]>>[R2X]-[C;!$(C[*;-1]);!$(C(=O)O);+0:1]-[R4X]",
                "description": "nucleophilic addition into electrophilic alpha carbon",
            },
            {
                "template": "[R1X].[C;!$(C[*;-1]);!$(C(=O)O);+0:1]=[C;+0:2]-[C&+0:3]=[R3X]>>[R2X]-[C;!$(C[*;-1]);!$(C(=O)O);+0:1]-[C;-1:2]-[C&+0:3]=[R3X]",
                "description": "nucleophilic addition into beta carbon1",
            },
            {
                "template": "[R2X]-[C;!$(C[*;-1]);!$(C(=O)O);+0:1]-[C;-1:2]-[C&+0:3]=[R3X]>>[R2X]-[C;!$(C[*;-1]);!$(C(=O)O);+0:1]-[C;+0:2]=[C&+0:3]-[R4X]",
                "description": "electron propagation",
            },
            # {
            #     "template": "[R3X]=[C;!$(C[*;-1]);!$(C(=O)O);+0:1][*:8][*;!a:2][*;!a:6][R1X]>>[R4X]-[C;!$(C[*;-1]);!$(C(=O)O);+0:1]1[*:8][*:2][*:6][R2X]1",
            #     "description": "5-mem nucleophilic addition into electrophilic alpha carbon",
            # },
            # {
            #     "template": "[R3X]=[C;!$(C[*;-1]);!$(C(=O)O);+0:1][*:9][*:8][*:2][*:6][R1X]>>[R4X]-[C;!$(C[*;-1]);!$(C(=O)O);+0:1]1[*:9][*:8][*:2][*:6][R2X]1",
            #     "description": "6-mem nucleophilic addition into electrophilic alpha carbon",
            # },
            # {
            #     "template": "[R3X]=[C;!$(C[*;-1]);!$(C(=O)O);+0:1][*:9][*:8][*:2][*:6][*:3][R1X]>>[R4X]-[C;!$(C[*;-1]);!$(C(=O)O);+0:1]1[*:9][*:8][*:2][*:6][*:3][R2X]1",
            #     "description": "7-mem nucleophilic addition into electrophilic alpha carbon",
            # },
            # {
            #     "template": "[R3X]=[C;!$(C[*;-1]);!$(C(=O)O);+0:1][*:9][*:8][*:2][*:6][*:3][*:7][R1X]>>[R4X]-[C;!$(C[*;-1]);!$(C(=O)O);+0:1]1[*:9][*:8][*:2][*:6][*:3][*:7][R2X]1",
            #     "description": "8-mem nucleophilic addition into electrophilic alpha carbon",
            # },
        ],
    },
    "Enamine Additions": {
        "scope": [
            {
                "name": "amine/thiol into iminium/sulfonium/oxonium",
                "[R1X]": "[N,S;+0:5]",
                "[R2X]": "[N,S;+1:5]",
                "[R3X]": "[N,S,O;+1:4]",
                "[R4X]": "[N,S,O;+0:4]",
            },
            {
                "name": "carbanion/hydroxide/alkoxide/thiolate into iminium/sulfonium/oxonium",
                "[R1X]": "[C,O,S;!$(O[N+]);-1:5]",
                "[R2X]": "[C,O,S;!$(O[N+]);+0:5]",
                "[R3X]": "[N,S,O;+1:4]",
                "[R4X]": "[N,S,O;+0:4]",
            },
            {
                "name": "amine/thiol into carbonyl/thial",
                "[R1X]": "[N,S;+0:5]",
                "[R2X]": "[N,S;+1:5]",
                "[R3X]": "[O,S;+0:4]",
                "[R4X]": "[O,S;-1:4]",
            },
            {
                "name": "carbanion/hydroxide/alkoxide/thiolate into carbonyl/thial",
                "[R1X]": "[C,O,S;!$(O[N+]);-1:5]",
                "[R2X]": "[C,O,S;!$(O[N+]);+0:5]",
                "[R3X]": "[O,S;+0:4]",
                "[R4X]": "[O,S;-1:4]",
            },
        ],
        "mechanism": [
            {
                "template": "[R1X]-[C;H0;+0:6]=[C;H1;+0:7]>>[R2X]=[C;H0;+0:6]-[C;H1;-1:7]",
                "description": "enamine addition into beta carbon",
            },
            {
                "template": "[R2X]=[C;H0;+0:6]-[C;H1;-1:7].[C;+0:1]=[C;+0:2]-[C&+0:3]=[R3X]>>[R2X]=[C;H0;+0:6]-[C;H1;+0:7]-[C;+0:1]-[C;-1:2]-[C&+0:3]=[R3X]",
                "description": "enamine propagations 1",
            },
            {
                "template": "[R2X]=[C;H0;+0:6]-[C;H1;+0:7]-[C;+0:1]-[C;-1:2]-[C&+0:3]=[R3X]>>[R2X]=[C;H0;+0:6]-[C;H1;+0:7]-[C;+0:1]-[C;+0:2]=[C&+0:3]-[R4X]",
                "description": "enamine propagations 2",
            },
        ],
    },
    "SN1": {
        "scope": [
            {
                "name": "halide/o-sub + tertiary carbon",
                "[R1X]": "[Cl,Br,I,O$(OS);H0;+0:2]",
                "[R2X]": "[Cl,Br,I,O$(OS);H0;-1:2]",
                "[R3X]": "[C&H0&+0:1]",
                "[R4X]": "[C&H0&+1:1]",
            },
        ],
        "mechanism": [
            {"template": "[R3X]-[R1X]>>[R4X].[R2X]", "description": "SN1 elimination"}
        ],
    },
    "SN2": {
        "scope": [
            {
                "name": "halide/O-sub + enamine",
                "[R1X]": "[Br,Cl,I,O;H0;+0:2]",
                "[R2X]": "[Br,Cl,I,O;H0;-1:2]",
                "[R3X]": "[N;+0:5]-[C;+0:6]=[C;+0:7]",
                "[R4X]": "[N;+1:5]=[C;+0:6]-[C;+0:7]",
            },
            {
                "name": "halide/O-sub + amine",
                "[R1X]": "[Br,Cl,I,O;H0;+0:2]",
                "[R2X]": "[Br,Cl,I,O;H0;-1:2]",
                "[R3X]": "[N,P&X3;+0:5]",
                "[R4X]": "[N,P&X3;+1:5]",
            },
            {
                "name": "halide + carbanion",
                "[R1X]": "[Br,Cl,I,O;H0;+0:2]",
                "[R2X]": "[Br,Cl,I,O;H0;-1:2]",
                "[R3X]": "[C,O,S;-1:5]",
                "[R4X]": "[C,O,S;+0:5]",
            },
        ],
        "mechanism": [
            {
                "template": "[R3X].[C;H3,H2,H1;+0:1]-[R1X]>>[R4X]-[C;H3,H2,H1;+0:1].[R2X]",
                "description": "SN2",
            },
            {
                "template": "[R3X].[c;+0:1]-[R1X]>>[R4X]-[c;+0:1].[R2X]",
                "description": "sNar",
            },
            # {
            #     "template": "[R3X][*:10][*:11][*:12][C;H2,H1;$(Caa[!a][C,O,S,N,P&H1]),$(Ca[!a][!a][C,O,S,N,P&H1]);+0:1]-[R1X]>>[R4X]([*:10][*:11][*:12]-1)-[C,H2,H1;$(Caa[!a][C,O,S,N,P&H1]),$(Ca[!a][!a][C,O,S,N,P&H1]);+0:1]1.[R2X]",
            #     "description": "5 mem intramolecular SN2",
            # },
            # {
            #     "template": "[R3X][*:10][*:11][*:12][*:13][C;H2,H1;+0:1]-[R1X]>>[R4X]([*:10][*:11][*:12][*:13]-1)-[C,H2,H1;+0:1]1.[R2X]",
            #     "description": "6 mem intramolecular SN2",
            # },
            # {
            #     "template": "[R3X][*:10][*:11][*:12][*:13][*:14][C;H2,H1;+0:1]-[R1X]>>[R4X]([*:10][*:11][*:12][*:13][*:14]-1)-[C,H2,H1;+0:1]1.[R2X]",
            #     "description": "7 mem intramolecular SN2",
            # },
            # {
            #     "template": "[R3X][*:10][*:11][*:12][*:13][*:14][*:15][C;H2,H1;+0:1]-[R1X]>>[R4X]([*:10][*:11][*:12][*:13][*:14][*:15]-1)-[C,H2,H1;+0:1]1.[R2X]",
            #     "description": "8 mem intramolecular SN2",
            # },
        ],
    },
    "Beta-elimination (intramolecular)": {
        "scope": [
            {
                "name": "primary or secondary beta carbon",
                "[R1X]": "[C;H2,H1;+0:1]",
                "[R2X]": "[C;H1,H0;+0:1]",
                "[R3X]": "[C&+0:2]",
            },
        ],
        "mechanism": [
            {
                "template": "[R1X]-[R3X]-[O&H0&+0:3]-[C&H0&+0:4]=[O&H0&+0:5]>>([R2X]=[R3X].[O&H0&+0:3]=[C&H0&+0:4]-[O&H1&+0:5])",
                "description": "intramolecular beta-elimination",
            }
        ],
    },
    "Arbuzov reaction": {
        "scope": [
            {
                "name": "carbon phosphite, primary or secondary halide",
                "[R1X]": "[C;+0:7]",
                "[R2X]": "[Cl,Br,I&H0&+0:6]",
                "[R3X]": "[Cl,Br,I&H0&-1:6]",
                "[R4X]": "[C;H2,H1;+0:5]",
            },
        ],
        "mechanism": [
            {
                "template": "[P&H0&+0:1](-[O&H0&+0:2][R1X])(-[O&H0&+0:3])-[O&H0&+0:4].[R4X]-[R2X]>>[P&H0&+:1](-[O&H0&+0:2][R1X])(-[O&H0&+0:3])(-[O&H0&+0:4])-[R4X].[R3X]",
                "description": "SN2 attack of alkyl halide",
            },
            {
                "template": "[P&H0&+:1](-[O&H0&+0:2][R1X])(-[O&H0&+0:3])(-[O&H0&+0:4])-[R4X].[R3X]>>[R2X][R1X].[P&H0&+0:1](=[O&H0&+0:2])(-[O&H0&+0:3])(-[O&H0&+0:4])-[R4X]",
                "description": "formation of phosphonate",
            },
        ],
    },
    "Sulfonium formation": {
        "scope": [
            {"name": "sulfonium", "[R1X]": "[S;+0:8]", "[R2X]": "[S;+1:8]"},
            {"name": "oxonium", "[R1X]": "[O;+0:8]", "[R2X]": "[O;+1:8]"},
        ],
        "mechanism": [
            {
                "template": "[R1X]-[C;!$(C(=O)O);+0:9](-[O;H2;+1:10])>>[R2X]=[C;!$(C(=O)O);+0:9].[O;H2;+0:10]",
                "description": "-ium formation via water extrusion",
            },
        ],
    },
    "iminium formation": {
        "scope": [
            {"name": "sulfonium", "[R1X]": "[S;+0:8]", "[R2X]": "[S;+1:8]"},
            {"name": "oxonium", "[R1X]": "[O;+0:8]", "[R2X]": "[O;+1:8]"},
            {"name": "iminium", "[R1X]": "[N;+0:8]", "[R2X]": "[N;+1:8]"},
        ],
        "mechanism": [
            {
                "template": "[S;+0:2][C;H2:3][R1X]>>[S;-1:2].[C;H2:3]=[R2X]",
                "description": "iminium formation via methylsulfur extrusion",
            },

        ],
    },
    "Ene Reaction": {
        "scope": [
            {
                "name": "ene",
                "[R1X]": "[C;H2,H1,H0;+0:1]",
                "[R2X]": "[C;H3,H2,H1;+0:1]",
            },
            {
                "name": "imine",
                "[R1X]": "[N;H0;+0:1]",
                "[R2X]": "[N;H1;+0:1]",
            },
            {
                "name": "ketone",
                "[R1X]": "[O;H0;!$(O=CN);+0:1]",
                "[R2X]": "[O;H1;!$(O=CN);+0:1]",
            },
        ],
        "mechanism": [
            {
                "template": "[R1X]=[C;+0:2].[C;+0:3]=[C;+0:4]-[C;H3,H2,H1;+0:5]>>[R2X]-[C;+0:2]-[C;+0:3]-[C;+0:4]=[C;H2,H1,H0;+0:5]",
                "description": "1p alder ene reaction",
            },
            # {
            #     "template": "[R1X]=[C;+0:2][*:6][*:7][*:8][C;+0:3]=[C;+0:4]-[C;H3,H2,H1;+0:5]>>[R2X]-[C;+0:2]1[*:6][*:7][*:8][C;+0:3]-1-[C;+0:4]=[C;H2,H1,H0;+0:5]",
            #     "description": "intramolecular 5-mem alder ene reaction",
            # },
            # {
            #     "template": "[R1X]=[C;+0:2][*:6][*:7][*:8][*:9][C;+0:3]=[C;+0:4]-[C;H3,H2,H1;+0:5]>>[R2X]-[C;+0:2]1[*:6][*:7][*:8][*:9][C;+0:3]-1-[C;+0:4]=[C;H2,H1,H0;+0:5]",
            #     "description": "intramolecular 6-mem alder ene reaction",
            # },
            # {
            #     "template": "[R1X]=[C;+0:2][*:6][*:7][*:8][*:10][*:9][C;+0:3]=[C;+0:4]-[C;H3,H2,H1;+0:5]>>[R2X]-[C;+0:2]1[*:6][*:7][*:8][*:10][*:9][C;+0:3]-1-[C;+0:4]=[C;H2,H1,H0;+0:5]",
            #     "description": "intramolecular 7-mem alder ene reaction",
            # },
            # {
            #     "template": "[R1X]=[C;+0:2][*:6][*:7][*:8][*:10][*:11][*:9][C;+0:3]=[C;+0:4]-[C;H3,H2,H1;+0:5]>>[R2X]-[C;+0:2]1[*:6][*:7][*:8][*:10][*:11][*:9][C;+0:3]-1-[C;+0:4]=[C;H2,H1,H0;+0:5]",
            #     "description": "intramolecular 8-mem alder ene reaction",
            # },
        ],
    },
    "Conia-ene Reaction": {
        "scope": [
            {
                "name": "enol 2/1p vinyl 2/1/0p",
                "[R1X]": "[C;H2,H1;+0:3]",
                "[R2X]": "[C;H1,H0;+0:3]",
                "[R3X]": "[C;H2,H1,H0;+0:8]",
                "[R4X]": "[C;H3,H2,H1;+0:8]",
            },
        ],
        "mechanism": [
            {
                "template": "[C;H1,H0;!$(C(=O)[O&H1,N,Cl,Br,I,F]);+0:2](-[R1X])=[O&H0&+0:9]>>[C;H1,H0;!$(C(=O)[O&H1,N,Cl,Br,I,F]);+0:2](=[R2X])-[O&H1&+0:9]",
                "description": "conia ene tautomerization",
            },
            {
                "template": "[C;H1,H0;!$(C(=O)[O&H1,N,Cl,Br,I,F]);+0:2](=[R2X])-[O&H1&+0:9].[C&+0:6]-[C&+0:7]=[R3X]>>[C;H1,H0;!$(C(=O)[O&H1,N,Cl,Br,I,F]);+0:2](-[R2X]-[C&+0:7](-[C&+0:6])-[R4X])=[O&H1&+1:9].",
                "description": "conia ene cyclization",
            },
            # {
            #     "template": "[C;H0,H1;!$(C(=O)[O&H1,N,Cl,Br,I,F]);+0:2](=[R2X]-[C&+0:4]-[C&+0:5]-[C&+0:6]-[C&+0:7]=[R3X])-[O&H1&+0:9]>>[C;H0,H1;!$(C(=O)[O&H1,N,Cl,Br,I,F]);+0:2](-[R2X]1-[C&+0:4]-[C&+0:5]-[C&+0:6]-[C&+0:7]-1-[R4X])=[O&H0&+0:9]",
            #     "description": "5mem conia ene cyclization",
            # },
            # {
            #     "template": "[C;H0,H1;!$(C(=O)[O&H1,N,Cl,Br,I,F]);+0:2](=[R2X]-[C&+0:4]-[C&+0:5]-[C&+0:12]-[C&+0:6]-[C&+0:7]=[R3X])-[O&H1&+0:9]>>[C;H0,H1;!$(C(=O)[O&H1,N,Cl,Br,I,F]);+0:2](-[R2X]1-[C&+0:4]-[C&+0:5]-[C&+0:12]-[C&+0:6]-[C&+0:7]-1-[R4X])=[O&H0&+0:9]",
            #     "description": "6mem conia ene cyclization",
            # },
            # {
            #     "template": "[C;H0,H1;!$(C(=O)[O&H1,N,Cl,Br,I,F]);+0:2](=[R2X]-[C&+0:4]-[C&+0:5]-[C&+0:10]-[C&+0:12]-[C&+0:6]-[C&+0:7]=[R3X])-[O&H1&+0:9]>>[C;H0,H1;!$(C(=O)[O&H1,N,Cl,Br,I,F]);+0:2](-[R2X]1-[C&+0:4]-[C&+0:5]-[C&+0:10]-[C&+0:12]-[C&+0:6]-[C&+0:7]-1-[R4X])=[O&H0&+0:9]",
            #     "description": "7mem conia ene cyclization",
            # },
            # {
            #     "template": "[C;H0,H1;!$(C(=O)[O&H1,N,Cl,Br,I,F]);+0:2](=[R2X]-[C&+0:4]-[C&+0:5]-[C&+0:11]-[C&+0:10]-[C&+0:12]-[C&+0:6]-[C&+0:7]=[R3X])-[O&H1&+0:9]>>[C;H0,H1;!$(C(=O)[O&H1,N,Cl,Br,I,F]);+0:2](-[R2X]1-[C&+0:4]-[C&+0:5]-[C&+0:11]-[C&+0:10]-[C&+0:12]-[C&+0:6]-[C&+0:7]-1-[R4X])=[O&H0&+0:9]",
            #     "description": "8mem conia ene cyclization",
            # },
        ],
    },
    "Prins Reaction": {
        "scope": [
            {
                "name": "prins aldehyde/acetone",
                "[R1X]": "[C;!$(C(=O)N);+0:1]",
                "[R2X]": "[C;!$(C(=O)N);+1:1]",
            },
        ],
        "mechanism": [
            {
                "template": "[R1X]=[O;H0;+0:2]>>[R1X]=[O;H1;+1:2]",
                "description": "protonation of aldehyde/acetone",
            },
            {
                "template": "[R1X]=[O;H1;+1:2]>>[R2X]-[O;H1;+0:2]",
                "description": "tautomerization of aldehyde/acetone",
            },
            {
                "template": "[R2X]-[O;H1;+0:2].[C;+0:3]=[C;H1,H0;+0:4]>>[R1X](-[C;+0:3]-[C;H1,H0;+1:4])-[O;H1;+0:2]",
                "description": "ene addition into carbocation",
            },
            {
                "template": "[R1X](-[C;+0:3]-[C;H1,H0;+1:4])-[O;H1;+0:2].[O;H2;+0:5]>>[R1X](-[C;+0:3]-[C;H1,H0;+0:4]-[O;H1;+0:5])-[O;H1;+0:2]",
                "description": "prins quench with water",
            },
            {
                "template": "[R1X](-[C;H2,H1;+0:3]-[C;H1,H0;+1:4])-[O;H1;+0:2]>>[R1X](-[C;H1,H0;+0:3]=[C;H1,H0;+0:4])-[O;H1;+0:2]",
                "description": "prins quench with base",
            },
            {
                "template": "[R1X](-[C;+0:3]-[C;H1,H0;+1:4])-[O;H1;+0:2].[C;+0:5]=[O;H0;+0:6]>>[R1X](-[C;+0:3]-[C;H1,H0;+0:4]([O;H0;+1:6]=[C;+0:5]))-[O;H1;+0:2]",
                "description": "prins carbonyl addition",
            },
            {
                "template": "[R1X](-[C;+0:3]-[C;H1,H0;+0:4]([O;H0;+1:6]=[C;+0:5]))-[O;H1;+0:2]>>[R1X](-[C;+0:3]-[C;H1,H0;+0:4]([O;H0;+0:6]-[C;+1:5]))-[O;H1;+0:2]",
                "description": "prins carbonyl addition tautomer",
            },
            {
                "template": "[R1X](-[C;+0:3]-[C;H1,H0;+0:4]([O;H0;+0:6]-[C;+1:5]))-[O;H1;+0:2]>>[R1X]1-[C;+0:3]-[C;H1,H0;+0:4]-[O&H0&+0:6]-[C;+0:5]-[O&H1&+:2]-1",
                "description": "prins carbonyl cyclization to gemdiol",
            },
            {
                "template": "[R1X]1-[C;+0:3]-[C;H1,H0;+0:4]-[O&H0&+0:6]-[C;+0:5]-[O&H1&+:2]-1>>[R1X]1-[C;+0:3]-[C;H1,H0;+0:4]-[O&H0&+0:6]-[C;+0:5]-[O&H0&+0:2]-1",
                "description": "deprotonation of gemdiol",
            },
        ],
    },
    "Pinacol Rearrangement": {
        "scope": [
            {"name": "carbon migration", "[R1X]": "[C;+0:5]"},
        ],
        "mechanism": [
            {
                "template": "[O;H1;+0:1]-[C;+0:2]-[C;+0:3](-[R1X])-[O;H1;+0:4]>>[O;H2;+1:1]-[C;+0:2]-[C;+0:3](-[R1X])-[O;H1;+0:4]",
                "description": "protonation of pinacol",
            },
            {
                "template": "[O;H2;+1:1]-[C;+0:2]-[C;+0:3](-[R1X])-[O;H1;+0:4]>>[O;H2;+0:1].[C;+1:2]-[C;+0:3](-[R1X])-[O;H1;+0:4]",
                "description": "loss of water from pinacol",
            },
            {
                "template": "[C;+1:2]-[C;+0:3](-[R1X])-[O;H1;+0:4]>>[C;+0:2](-[R1X])-[C;+1:3]-[O;H1;+0:4]",
                "description": "carbocation migration",
            },
            {
                "template": "[C;+0:2](-[R1X])-[C;+1:3]-[O;H1;+0:4]>>[C;+0:2](-[R1X])-[C;+0:3]=[O;H0;+0:4]",
                "description": "deprotonation of alcohol to quench carbocation",
            },
        ],
    },
    "Benzilic Acid Rearrangement": {
        "scope": [
            {
                "name": "hydroxyl nucleophile",
                "[R1X]": "[O;!$(O[N+]);H0;-1:5]",
                "[R2X]": "[O;!$(O[N+]);H0;+0:5]",
            },
            {
                "name": "amide anion nucleophile",
                "[R1X]": "[N;H0,H1,H2;-1:5]",
                "[R2X]": "[N;H0,H1,H2;+0:5]",
            },
        ],
        "mechanism": [
            {
                "template": "[O;H0;+0:1]=[C;H0;+0:2]([#6;+0:6])-[C;H0;+0:3]=[O;H0;+0:4].[R1X]>>[O;H0;-1:1]-[C;H0;+0:2]([#6;+0:6])(-[R2X])-[C;H0;+0:3]=[O;H0;+0:4]",
                "description": "nucleophilic attack of diketone",
            },
            {
                "template": "[O;H0;-1:1]-[C;H0;+0:2]([#6;+0:6])(-[R2X])-[C;H0;+0:3]=[O;H0;+0:4]>>[O;H0;+0:1]=[C;H0;+0:2](-[R2X])-[C;H0;+0:3]([#6;+0:6])-[O;H0;-1:4]",
                "description": "alkyl migration",
            },
            {
                "template": "[O;H0;+0:1]=[C;H0;+0:2](-[R2X])-[C;H0;+0:3]([#6;+0:6])-[O;H0;-1:4]>>[O;H0;+0:1]=[C;H0;+0:2](-[R2X])-[C;H0;+0:3]([#6;+0:6])-[O;H1;+0:4]",
                "description": "protonation",
            },
        ],
    },
    "Cannizzaro Reaction": {
        "scope": [
            {
                "name": "hydroxyl/alkoxide nucleophile",
                "[R1X]": "[O;!$(O[N+]);H1,H0;-1:5]",
                "[R2X]": "[O;!$(O[N+]);H1,H0;+0:5]",
            },
        ],
        "mechanism": [
            {
                "template": "[O;H0;+0:1]=[C;H1;+0:2].[R1X]>>[O;H0;-1:1]-[C;H1;+0:2]-[R2X]",
                "description": "nucleophilic attack of aldehyde",
            },
            {
                "template": "[O;H0;-1:1]-[C;H1;+0:2]-[R2X].[O;H0;+0:3]=[C;H1;+0:4]>>[O;H0;+0:1]=[C;H0;+0:2]-[R2X].[O;H0;-1:3]-[C;H2;+0:4]",
                "description": "hydride transfer",
            },
        ],
    },
    "Claisen Rearrangement": {
        "scope": [
            {
                "name": "vinyl and ether",
                "[R2X]": "[O&H0&+0:1]",
                "[R3X]": "[O&H1&+0:1]",
            },
            {
                "name": "vinyl and azaclaisen p0 or p1",
                "[R2X]": "[N;H0,H1;+0:1]",
                "[R3X]": "[N;H1,H2;+0:1]",
            },
        ],
        "mechanism": [
            {
                "template": "[R2X](-[C&+0:2]-[C&H1&+0:3]=[C;H2,H1;+0:4])-[C&+0:5]=[C&+0:6]>>[R2X]=[C&+0:5]-[C&+0:6]-[C;H2,H1;+0:4]-[C&H1&+0:3]=[C&+0:2]",
                "description": "p1/p2 'vinyl' claisen rearrangement",
            },
            {
                "template": "[R2X](-[C&+0:2]-[C&H1&+0:3]=[C;H2,H1;+0:4])-[c;H0;+0:5][c;H1;+0:6]>>[R3X]-[c;H0;+0:5][c;H0;+0:6]-[C;H2,H1;+0:4]-[C&H1&+0:3]=[C&+0:2]",
                "description": "p1/p2 aromatic claisen rearrangement",
            },
        ],
    },
    "Eschenmoser-Claisen Rearrangement": {
        "scope": [
            {
                "name": "p3 or p2",
                "[R1X]": "[C;H3,H2;+0:11]",
                "[R2X]": "[C;H2,H1;+0:11]",
            },
        ],
        "mechanism": [
            {
                "template": "[N&H0&+0:1]-[C&H0&+0:2](-[R1X])(-[O&H0&+0:3])-[O&H0&+0:4]>>([N&H0&+:1]=[C&H0&+0:2](-[R1X])-[O&H0&+0:4].[O&H0&-1:3])",
                "description": "iminium formation by loss of alkoxide",
            },
            {
                "template": "[N&H0&+:1]=[C&H0&+0:2](-[R1X])-[O&H0&+0:4].[O&H1&+0:5]-[C&+0:6](-[C&+0:7])-[C&+0:8]=[C&+0:9]-[C&+0:10]>>[N&H0&+0:1]-[C&H0&+0:2](-[R1X])(-[O&H0&+0:4])-[O&H1&+:5]-[C&+0:6](-[C&+0:7])-[C&+0:8]=[C&+0:9]-[C&+0:10]",
                "description": "addition of enol",
            },
            {
                "template": "[N&H0&+0:1]-[C&H0&+0:2](-[R1X])(-[O&H0&+0:4])-[O&H1&+:5]-[C&+0:6](-[C&+0:7])-[C&+0:8]=[C&+0:9]-[C&+0:10]>>[N&H0&+0:1]-[C&H0&+0:2](-[R1X])(-[O&H1&+1:4])-[O&H0&+0:5]-[C&+0:6](-[C&+0:7])-[C&+0:8]=[C&+0:9]-[C&+0:10]",
                "description": "proton shift",
            },
            {
                "template": "[N&H0&+0:1]-[C&H0&+0:2](-[R1X])(-[O&H1&+1:4])-[O&H0&+0:5]-[C&+0:6](-[C&+0:7])-[C&+0:8]=[C&+0:9]-[C&+0:10]>>([N&H0&+0:1]-[C&H0&+0:2](=[R2X])-[O&H0&+0:5]-[C&+0:6](-[C&+0:7])-[C&+0:8]=[C&+0:9]-[C&+0:10].[O&H1&+0:4])",
                "description": "loss of second alkoxide",
            },
            {
                "template": "[N&H0&+0:1]-[C&H0&+0:2](=[R2X])-[O&H0&+0:5]-[C&+0:6](-[C&+0:7])-[C&+0:8]=[C&+0:9]-[C&+0:10]>>[N&H0&+0:1]-[C&H0&+0:2](=[O&H0&+0:5])-[R2X]-[C&+0:9](-[C&+0:8]=[C&+0:6]-[C&+0:7])-[C&+0:10]",
                "description": "claisen rearrangement",
            },
        ],
    },
    "Johnson-Claisen Rearrangement": {
        "scope": [
            {
                "name": "methyl/alkyl",
                "[R1X]": "[C;H3,H2;+0:1]",
                "[R2X]": "[C;H2,H1;+0:1]",
            },
        ],
        "mechanism": [
            {
                "template": "[R1X]-[C&H0&+0:2](-[O&H0&+0:3])(-[O&H0&+0:4])-[O&H0&+0:5]>>[R1X]-[C&H0&+0:2](-[O&H1&+1:3])(-[O&H0&+0:4])-[O&H0&+0:5]",
                "description": "protonation of alkoxide",
            },
            {
                "template": "[R1X]-[C&H0&+0:2](-[O&H1&+1:3])(-[O&H0&+0:4])-[O&H0&+0:5]>>[R1X]-[C&H0&+0:2](=[O&H0&+1:4])-[O&H0&+0:5].[O&H1&+0:3]",
                "description": "loss of alkoxide",
            },
            {
                "template": "[R1X]-[C&H0&+0:2](=[O&H0&+:4])-[O&H0&+0:5].[C&H0&+0:7](-[C&H1&+0:8]=[C&H1&+0:9]-[C&+0:10])-[O&H1&+0:11]>>[R1X]-[C&H0&+0:2](-[O&H0&+0:4])(-[O&H0&+0:5])-[O&H1&+:11]-[C&H0&+0:7]-[C&H1&+0:8]=[C&H1&+0:9]-[C&+0:10]",
                "description": "addition of alcohol",
            },
            {
                "template": "[R1X]-[C&H0&+0:2](-[O&H0&+0:4])(-[O&H0&+0:5])-[O&H1&+:11]-[C&H0&+0:7]-[C&H1&+0:8]=[C&H1&+0:9]-[C&+0:10]>>[R1X]-[C&H0&+0:2](-[O&H0&+0:4])(-[O&H1&+:5])-[O&H0&+0:11]-[C&H0&+0:7]-[C&H1&+0:8]=[C&H1&+0:9]-[C&+0:10]",
                "description": "proton transfer",
            },
            {
                "template": "[R1X]-[C&H0&+0:2](-[O&H0&+0:4])(-[O&H1&+:5])-[O&H0&+0:11]-[C&H0&+0:7]-[C&H1&+0:8]=[C&H1&+0:9]-[C&+0:10]>>[R2X]=[C&H0&+0:2](-[O&H0&+0:4])-[O&H0&+0:11]-[C&H0&+0:7]-[C&H1&+0:8]=[C&H1&+0:9]-[C&+0:10].[O&H1&+0:5]",
                "description": "loss of second alkoxide",
            },
            {
                "template": "[R2X]=[C&H0&+0:2](-[O&H0&+0:4])-[O&H0&+0:11]-[C&H0&+0:7]-[C&H1&+0:8]=[C&H1&+0:9]-[C&+0:10]>>[R2X](-[C&H0&+0:2](-[O&H0&+0:4])=[O&H0&+0:11])-[C&H1&+0:9](-[C&H1&+0:8]=[C&H0&+0:7])-[C&+0:10]",
                "description": "claisen rearrangement",
            },
        ],
    },
    "Carroll Rearrangement": {
        "scope": [
            {"name": "ether", "[R1X]": "[O&H0&+0:4]"},
        ],
        "mechanism": [
            {
                "template": "[C&H2&+0:1]=[C&H1&+0:2]-[C&H2&+0:3]-[R1X]-[C&H0&+0:5](=[O&H0&+0:6])-[C&H2&+0:7]-[C&H0&+0:8]=[O&H0&+0:9]>>[C&H2&+0:1]=[C&H1&+0:2]-[C&H2&+0:3]-[R1X]-[C&H0&+0:5](-[O&H1&+0:6])=[C&H1&+0:7]-[C&H0&+0:8]=[O&H0&+0:9]",
                "description": "base mediated formation of enol",
            },
            {
                "template": "[C&H2&+0:1]=[C&H1&+0:2]-[C&H2&+0:3]-[R1X]-[C&H0&+0:5](-[O&H1&+0:6])=[C&H1&+0:7]-[C&H0&+0:8]=[O&H0&+0:9]>>[C&H2&+0:1](-[C&H1&+0:2]=[C&H2&+0:3])-[C&H1&+0:7](-[C&H0&+0:5](=[R1X])-[O&H1&+0:6])-[C&H0&+0:8]=[O&H0&+0:9]",
                "description": "rearrangement of enol",
            },
            {
                "template": "[C&H2&+0:1](-[C&H1&+0:2]=[C&H2&+0:3])-[C&H1&+0:7](-[C&H0&+0:5](=[R1X])-[O&H1&+0:6])-[C&H0&+0:8]=[O&H0&+0:9]>>[C&H2&+0:1](-[C&H1&+0:2]=[C&H2&+0:3])-[C&H2&+0:7]-[C&H0&+0:8]=[O&H0&+0:9].[R1X]=[C&H0&+0:5]=[O&H0&+0:6]",
                "description": "loss of CO2",
            },
        ],
    },
    "Wolff Rearrangement": {
        "scope": [
            {"name": "alkyl", "[R1X]": "[#6&+0:3]"},
        ],
        "mechanism": [
            {
                "template": "[O&H0&+0:1]=[C&H0&+0:2](-[R1X])-[C&H0&+0:4]=[N&H0&+:5]=[N&H0&-:6]>>[O&H0&-1:1]-[C&H0&+0:2](-[R1X])=[C&H0&+0:4]-[N&H0&+:5]#[N&H0&+0:6]",
                "description": "diazo tautomerization",
            },
            {
                "template": "[O&H0&-1:1]-[C&H0&+0:2](-[R1X])=[C&H0&+0:4]-[N&H0&+:5]#[N&H0&-:6]>>[O&H0&+0:1]=[C&H0&+0:2]=[C&H0&+0:4](-[R1X]).[N&H0&+0:5]#[N&H0&+0:6]",
                "description": "1,2 rearrangement to form ketene",
            },
        ],
    },
    "Westphalen–Lettre Rearrangement": {
        "scope": [
            {"name": "alkyl", "[R1X]": "[#6&+0:3]"},
        ],
        "mechanism": [
            {
                "template": "[C&+0:1]1-[C&H0&+0:2]2(-[C&H1&+0:3]-[C&+0:4]-[C&+0:7]-[C&+0:9]-[C&H0&+0:10]-2(-[C&+0:8]-[C&+0:6]-[C&+0:5]-1)-[O&H1&+0:12])-[#6&+0:11]>>[C&+0:1]1-[C&H0&+0:2]2(-[C&H1&+0:3]-[C&+0:4]-[C&+0:7]-[C&+0:9]-[C&H0&+:10]-2-[C&+0:8]-[C&+0:6]-[C&+0:5]-1)-[#6&+0:11].[O&H2&+0:12]",
                "description": "loss of water to form carbocation",
            },
            {
                "template": "[C&+0:1]1-[C&H0&+0:2]2(-[C&H1&+0:3]-[C&+0:4]-[C&+0:7]-[C&+0:9]-[C&H0&+:10]-2-[C&+0:8]-[C&+0:6]-[C&+0:5]-1)-[#6&+0:11]>>[C&+0:1]1-[C&H0&+0:2]2=[C&H0&+0:3]-[C&+0:4]-[C&+0:7]-[C&+0:9]-[C&H0&+:10]-2(-[C&+0:8]-[C&+0:6]-[C&+0:5]-1)-[#6&+0:11]",
                "description": "methyl rearrangement",
            },
        ],
    },
    "Robinson Annulation": {
        "scope": [
            {"name": "hydroxide", "[R1X]": "[O&H1&-:11]"},
        ],
        "mechanism": [
            {
                "template": "[C&+0:1]1-[C&H0&+0:2]2(-[C&+0:3]=[C&+0:4](-[C&+0:7]-[C&+0:9]-[C&+0:10]-2-[C&+0:8]-[C&+0:6]-[C&+0:5]-1)-[O&H0&-:12])-[O&H1&+0:11]>>[C&+0:1]1-[C&H0&+0:2]2=[C&+0:3]-[C&+0:4](-[C&+0:7]-[C&+0:9]-[C&+0:10]-2-[C&+0:8]-[C&+0:6]-[C&+0:5]-1)=[O&H0&+0:12].[R1X]",
                "description": "loss of water to finish annulation",
            },
        ],
    },
    "Bellus–Claisen Rearrangement": {
        "scope": [
            {
                "name": "ether",
                "[R2X]": "[O&H0&+0:1]",
                "[R3X]": "[O&H0&+:1]",
            },
            {
                "name": "thioether",
                "[R2X]": "[S&H0&+0:1]",
                "[R3X]": "[S&H0&+:1]",
            },
            {
                "name": "amine",
                "[R2X]": "[N&+0:1]",
                "[R3X]": "[N&+:1]",
            },
        ],
        "mechanism": [
            {
                "template": "[R2X].[C&H0&+0:8]=[C&H0&+0:9]=[O&H0&+0:10]>>[R3X]-[C&H0&+0:9](=[C&H0&+0:8])-[O&H0&-:10]",
                "description": "addition of ketene",
            },
            {
                "template": "[R3X]1(-[C&H1&+0:2](-[C&H0&+0:3]=[C&H2&+0:4])-[C&H2&+0:6]-[C&H2&+0:5]-1)-[C&H0&+0:9](=[C&H0&+0:8])-[O&H0&-:10]>>[C&H1&+0:2](-[C&H2&+0:6]-[C&H2&+0:5]-[R2X]1)=[C&H0&+0:3]-[C&H2&+0:4]-[C&H0&+0:8]-[C&H0&+0:9]-1=[O&H0&+0:10]",
                "description": "4 mem ring expansion",
            },
            {
                "template": "[R3X]1(-[C&H1&+0:2](-[C&H0&+0:3]=[C&H2&+0:4])-[C&H2&+0:6]-[C&H2&+0:7]-[C&H2&+0:5]-1)-[C&H0&+0:9](=[C&H0&+0:8])-[O&H0&-:10]>>[C&H1&+0:2](-[C&H2&+0:6]-[C&H2&+0:7]-[C&H2&+0:5]-[R2X]1)=[C&H0&+0:3]-[C&H2&+0:4]-[C&H0&+0:8]-[C&H0&+0:9]-1=[O&H0&+0:10]",
                "description": "5 mem ring expansion",
            },
            {
                "template": "[R3X]1(-[C&H1&+0:2](-[C&H0&+0:3]=[C&H2&+0:4])-[C&H2&+0:6]-[C&H2&+0:7]-[C&H2&+0:11]-[C&H2&+0:5]-1)-[C&H0&+0:9](=[C&H0&+0:8])-[O&H0&-:10]>>[C&H1&+0:2](-[C&H2&+0:6]-[C&H2&+0:7]-[C&H2&+0:11]-[C&H2&+0:5]-[R2X]1)=[C&H0&+0:3]-[C&H2&+0:4]-[C&H0&+0:8]-[C&H0&+0:9]-1=[O&H0&+0:10]",
                "description": "6 mem ring expansion",
            },
        ],
    },
    "Alpha-ketol Rearrangement": {
        "scope": [
            {"name": "tertiary carbon", "[R1X]": "[C;H0;+0:3]"},
        ],
        "mechanism": [
            {
                "template": "[O;H0;+0:1]=[C;H0;+0:2]-[R1X]([#6;+0:6])-[O;H1;+0:4]>>[O;H0;+0:1]=[C;H0;+0:2]-[R1X]([#6;+0:6])-[O;H0;-1:4]",
                "description": "deprotonation of ketol",
            },
            {
                "template": "[O;H0;+0:1]=[C;H0;+0:2]-[R1X]([#6;+0:6])-[O;H0;-1:4]>>[O;H1;+0:1]-[C;H0;+0:2]([#6;+0:6])-[R1X]=[O;H0;+0:4]",
                "description": "alkyl migration",
            },
        ],
    },
    "Aldose-ketose Transformation": {
        "scope": [
            {
                "name": "secondary carbon",
                "[R1X]": "[C;H1;+0:3]",
                "[R2X]": "[C;H0;+0:3]",
            },
        ],
        "mechanism": [
            {
                "template": "[O;H0;+0:1]=[C;H1;+0:2]-[R1X]-[O;H1;+0:4]>>[O;H0;-1:1]-[C;H1;+0:2]=[R2X]-[O;H1;+0:4]",
                "description": "tautomerization of aldose",
            },
            {
                "template": "[O;H0;-1:1]-[C;H1;+0:2]=[R2X]-[O;H1;+0:4]>>[O;H1;+0:1]-[C;H1;+0:2]=[R2X]-[O;H0;-1:4]",
                "description": "proton transfer",
            },
            {
                "template": "[O;H1;+0:1]-[C;H1;+0:2]=[R2X]-[O;H0;-1:4]>>[O;H1;+0:1]-[C;H2;+0:2]-[R2X]=[O;H0;+0:4]",
                "description": "tautomerization of ketose",
            },
        ],
    },
    "Cope Rearrangement": {
        "scope": [
            {"name": "cope", "[R1X]": "[#6;+0:7]"},
            {"name": "oxy-cope", "[R1X]": "[O;H1;+0:7]"},
        ],
        "mechanism": [
            {
                "template": "[C&+0:1](-[C&+0:2]=[C&+0:3])(-[C&+0:4]-[C&+0:5]=[C&+0:6])-[R1X]>>[C&+0:1](=[C&+0:2]-[C&+0:3]-[C&+0:6]-[C&+0:5]=[C&+0:4])-[R1X]",
                "description": "cope rearrangement",
            },
            {
                "template": "[C&+0:1](-[N;H1;+1:2]=[C&+0:3])(-[C&+0:4]-[C&+0:5]=[C&+0:6])-[R1X]>>[C&+0:1](=[N;H1;+1:2]-[C&+0:3]-[C&+0:6]-[C&+0:5]=[C&+0:4])-[R1X]",
                "description": "aza-cope rearrangement",
            },
        ],
    },
    "Aza-Cope Mannich": {
        "scope": [
            {"name": "cope", "[R1X]": "[#6;+0:7]"},
        ],
        "mechanism": [
            {
                "template": "[C&+0:1]=[N;+1:2]-[C&+0:3](-[R1X])-[C;H1;+0:6](-[O;H1;+0:8])-[C&+0:5]=[C&+0:4]>>[C&+0:1](-[N;+1:2]=[C&+0:3](-[R1X]))-[C&+0:4]-[C&+0:5]=[C;H1;+0:6]-[O&H1&+0:8]",
                "description": "aza-cope rearrangement",
            },
            {
                "template": "[C&+0:1](-[N;+1:2]=[C&+0:3](-[R1X]))-[C&+0:4]-[C&+0:5]=[C;H1;+0:6]-[O&H1&+0:8]>>[C&+0:1]1-[N;+0:2]-[C&+0:3](-[R1X])-[C&+0:5](-[C;H1;+0:6]=[O&H0&+0:8])-[C&+0:4]-1",
                "description": "Mannich cyclization",
            },
        ],
    },
    "Aza-Prins Pinacol": {
        "scope": [
            {"name": "cope", "[R1X]": "[#6;+0:7]"},
        ],
        "mechanism": [
            {
                "template": "[C&+0:1]=[N;+1:2]-[C&+0:3](-[R1X])-[C;H1;+0:6](-[O;H1;+0:8])-[C&+0:5]=[C&+0:4]>>[C&+0:1]1-[N;+0:2]-[C&+0:3](-[R1X])-[C;H1;+0:6](-[O&H1&+0:8])-[C&+:5]-[C&+0:4]-1",
                "description": "aza-prins cyclization",
            },
            {
                "template": "[C&+0:1]1-[N;+0:2]-[C&+0:3](-[R1X])-[C;H1;+0:6](-[O&H1&+0:8])-[C&+:5]-[C&+0:4]-1>>[C&+0:1]1-[N;+0:2]-[C&+0:3](-[R1X])-[C&+0:5](-[C;H1;+1:6]-[O&H1&+0:8])-[C&+0:4]-1",
                "description": "pinacol rearrangement",
            },
            {
                "template": "[C&+0:1]1-[N;+0:2]-[C&+0:3](-[R1X])-[C&+0:5](-[C;H1;+1:6]-[O&H1&+0:8])-[C&+0:4]-1>>[C&+0:1]1-[N;+0:2]-[C&+0:3](-[R1X])-[C&+0:5](-[C;H1;+0:6]=[O&H0&+0:8])-[C&+0:4]-1",
                "description": "aza-prins pinacol quench",
            },
        ],
    },
    "Ketene formation": {
        "scope": [
            {
                "name": "halide",
                "[R1X]": "[Cl,Br,I;H0;+0:1]",
                "[R2X]": "[Cl,Br,I;H0;-1:1]",
            },
        ],
        "mechanism": [
            {
                "template": "[C;H1;+0:8]-[C;H0;+0:9](=[O;H0;+0:10])[R1X]>>[C;H0;+0:8]=[C;H0;+0:9](=[O;H0;+0:10]).[R2X]",
                "description": "elimination into ketene",
            },
        ],
    },
    "Ketene additions": {
        "scope": [
            {
                "name": "water/alcohol (forms ester/acid)",
                "[R1X]": "[O;H2,H1;+0:1]",
                "[R2X]": "[O;H1,H0;+0:1]",
            },
            {
                "name": "hydroxide/alkoxide (forms ester/acid)",
                "[R1X]": "[O;H1,H0;-1:1]",
                "[R2X]": "[O;H1,H0;+0:1]",
            },
            {
                "name": "amine (forms amide)",
                "[R1X]": "[N;H3,H2,H1;+0:1]",
                "[R2X]": "[N;H2,H1,H0;+0:1]",
            },
            {
                "name": "thiol (forms thioester)",
                "[R1X]": "[S;H1;+0:1]",
                "[R2X]": "[S;H0;+0:1]",
            },
        ],
        "mechanism": [
            {
                "template": "[C;H0;+0:8]=[C;H0;+0:9](=[O;H0;+0:10]).[R1X]>>[C;H1;+0:8][C;H0;+0:9]([R2X])(=[O;H0;+0:10])",
                "description": "elimination into ketene",
            },
        ],
    },
    "Cycloaddition (2+2)": {
        "scope": [
            {"name": "imine", "[R1X]": "[N&+0:1]", "[R2X]": "[N&+0:1]"},
            {"name": "alkene", "[R1X]": "[C;+0:1]", "[R2X]": "[C;+0:1]"},
        ],
        "mechanism": [
            {
                "template": "[C;H0;+0:8]=[C;H0;+0:9](=[O;H0;+0:10]).[C;+0:11]=[R1X]>>[C&H0&+0:8]1-[C&H0&+0:9](=[O&H0&+0:10])-[R2X]-[C&+0:11]-1",
                "description": "2+2 cycloaddition",
            },
        ],
    },
    # "SEAr": {
    #     "scope": [
    #         {"name": "carbocation", "[R1X]": "[C;H2,H1,H0;!$(C=[OH+]);+1:1]", "[R2X]": "[C;H2,H1,H0;+0:1]"},
    #     ],
    #     "mechanism": [
    #         {
    #             "template": "[c;H1;+0:8]:[c;H1;+0:9].[R1X]>>[c;H1;+1:8]:[c;H1;+0:9]-[R2X]",
    #             "description": "carbocation addition into aromatic ring",
    #         },
    #         {
    #             "template": "[c;H1;+1:8]:[c;H1;+0:9]-[R2X]>>[c;H1;+0:8]:[c;H0;+0:9]-[R2X]",
    #             "description": "deprotonation of aryl carbocation",
    #         },
    #     ],
    # },
    # "SNAr": {
    #     "scope": [
    #         {
    #             "name": "hydroxyl/thiol",
    #             "[R1X]": "[#6;+0:4]",
    #             "[R2X]": "[C,O,S;!$(O[N+]);-1:5]",
    #             "[R3X]": "[C,O,S;!$(O[N+]);+0:5]",
    #         },
    #         {
    #             "name": "amine",
    #             "[R1X]": "[#6;+0:4]",
    #             "[R2X]": "[N,P&H1;+0:5]",
    #             "[R3X]": "[N,P&H1;+1:5]",
    #         },
    #     ],
    #     "mechanism": [
    #         {
    #             "template": "[c;+0:1]:[c;+0:2]-[Cl,Br,I;H0;+0:3].[R1X]-[R2X]>>[c;-1:1]:[c;+0:2]([R3X]-[R1X])-[Cl,Br,I;H0;+0:3]",
    #             "description": "Addition to aromatic ring",
    #         },
    #         {
    #             "template": "[c;-1:1]:[c;+0:2]([R3X]-[R1X])-[Cl,Br,I;H0;+0:3]>>[c;+0:1]:[c;+0:2]([R3X]-[R1X]).[Cl,Br,I;H0;-1:3]",
    #             "description": "Halide leaving",
    #         },
    #         {
    #             "template": "[c;+0:1]:[c;+0:2](-[Cl,Br,I;H0;+0:3])[c;+0:6]-[*:8]-[R1X]-[R2X]>>[c&-1:1]:[c&+0:2]1(-[Cl,Br,I;H0;+0:3])-[c&+0:6]-[*:8]-[R1X]-[R3X]-1",
    #             "description": "5 mem addition to aromatic ring",
    #         },
    #         {
    #             "template": "[c&-1:1]:[c&+0:2]1(-[Cl,Br,I;H0;+0:3])-[c&+0:6]-[*:8]-[R1X]-[R3X]-1>>[c&+0:1]:[c&+0:2]1-[c&+0:6]-[*:8]-[R1X]-[R3X]-1.[Cl,Br,I;H0;-1:3]",
    #             "description": "5 mem halide leaving",
    #         },
    #     ],
    # },
    "Deprotonations": {
        "scope": [
            {
                "name": "carboxylic acid",
                "[R1X]": "[O;H1;$(O([C;H0;+0]=[O;H0;+0]));+0:4]",
                "[R2X]": "[O;H0;$(O([C;H0;+0]=[O;H0;+0]));-1:4]",
            },
            {
                "name": "protonated carbonyl",
                "[R1X]": "[O;H1;+1:4]",
                "[R2X]": "[O;H0;+0:4]",
            },
            {"name": "alcohol", "[R1X]": "[O;H1;+0:4]", "[R2X]": "[O;H0;-1:4]"},
            {"name": "thiol", "[R1X]": "[S;H1;+0:4]", "[R2X]": "[S;H0;-1:4]"},
            {"name": "phosphonic acid", "[R1X]": "[P;H1;+1:4]", "[R2X]": "[P;H0;+0:4]"},
            {
                "name": "amines",
                "[R1X]": "[#7;H3,H2,H1;+1:4]",
                "[R2X]": "[#7;H2,H1,H0;+0:4]",
            },
        ],
        "mechanism": [
            {
                "template": "[R1X]>>[R2X]",
                "description": "Deprotonation",
            },
        ],
    },
    "Hydantoin Formation": {
        "scope": [
            {"name": "urea", "[R1X]": "[O;H0;+0:1]", "[R2X]": "[O;H1;+0:1]"},
            {"name": "thiourea", "[R1X]": "[S;H0;+0:1]", "[R2X]": "[S;H1;+0:1]"},
        ],
        "mechanism": [
            {
                "template": "[N;H2;+0:2]-[C;H0;+0:3](-[N;H2;+0:4])=[R1X]>>[N;H1;+0:2]=[C;H0;+0:3](-[N;H2;+0:4])-[R2X]",
                "description": "conjugation of urea",
            },
            {
                "template": "[N;H2;+0:5][C;H1;+0:6][C;H0;+0:7](-[O;H1;+0:8])=[O;H0;+0:9].[N;H1;+0:2]=[C;H0;+0:3](-[N;H2;+0:4])-[R2X]>>[N;H2;+0:5][C;H1;+0:6][C;H0;+0:7](-[O;H1;+0:8])(-[O;H1;+0:9])-[N;H1;+0:2]-[C;H0;+0:3](-[N;H2;+0:4])=[R1X]",
                "description": "addition into alpha amino acid",
            },
            {
                "template": "[N;H2;+0:5][C;H1;+0:6][C;H0;+0:7](-[O;H1;+0:8])(-[O;H1;+0:9])-[N;H1;+0:2]-[C;H0;+0:3](-[N;H2;+0:4])=[R1X]>>[N;H2;+0:5][C;H1;+0:6][C;H0;+0:7](=[O;H0;+0:8])-[N;H1;+0:2]-[C;H0;+0:3](-[N;H2;+0:4])=[R1X].[O;H2;+0:9]",
                "description": "condensation of substituted alpha amino acid",
            },
            {
                "template": "[N;H2;+0:5][C;H1;+0:6][C;H0;+0:7](=[O;H0;+0:8])-[N;H1;+0:2]-[C;H0;+0:3](-[N;H2;+0:4])=[R1X]>>[N&H1&+0:5]1-[C&H1&+0:6]-[C&H0&+0:7](=[O&H0&+0:8])-[N&H1&+0:2]-[C&H0&+0:3]-1(-[N&H2&+0:4])-[R2X]",
                "description": "cyclization of peptide",
            },
            {
                "template": "[N&H1&+0:5]1-[C&H1&+0:6]-[C&H0&+0:7](=[O&H0&+0:8])-[N&H1&+0:2]-[C&H0&+0:3]-1(-[N&H2&+0:4])-[R2X]>>[N&H1&+0:5]1-[C&H1&+0:6]-[C&H0&+0:7](=[O&H0&+0:8])-[N&H1&+0:2]-[C&H0&+0:3]-1=[R1X].[N&H3&+0:4]",
                "description": "loss of ammonia",
            },
        ],
    },
    "Protonations": {
        "scope": [
            {
                "name": "hydroxyl",
                "[R1X]": "[O;H0;!$(OC[O;H0;+0]);!$(O[*;+]);-1:4]",
                "[R2X]": "[O;H1;!$(OC[O;H0;+0]);!$(O[*;+]);+0:4]",
            },
            {
                "name": "alcohol",
                "[R1X]": "[O;H1;+0:4]",
                "[R2X]": "[O;H2;+1:4]",
            },
            {
                "name": "carbonyl",
                "[R1X]": "[O;H0;$(O=C);+0:4]",
                "[R2X]": "[O;H1;$(O=C);+1:4]",
            },
            {"name": "thiol", "[R1X]": "[S;H0;-1:4]", "[R2X]": "[S;H1;+0:4]"},
            {
                "name": "thiomethanol",
                "[R1X]": "[S;+0:4][C:5][O;H1;+0:6]",
                "[R2X]": "[S;+0:4][C:5][O;H2;+1:6]",
            },
            {"name": "phosphonic acid", "[R1X]": "[P;H0;-1:4]", "[R2X]": "[P;H1;+0:4]"},
            {
                "name": "carbanion",
                "[R1X]": "[C;H2,H1,H0;-1:4]",
                "[R2X]": "[C;H3,H2,H1;+0:4]",
            },
            {
                "name": "amide anion",
                "[R1X]": "[N;H2,H1,H0;-1:4]",
                "[R2X]": "[N;H3,H2,H1;+0:4]",
            },
        ],
        "mechanism": [
            {
                "template": "[R1X]>>[R2X]",
                "description": "Protonation",
            },
        ],
    },
    "Protection (N-Boc)": {
        "scope": [
            {
                "name": "amines",
                "[R1X]": "[N;H2,H1;+0:2]",
                "[R2X]": "[N;H2,H1;+1:2]",
                "[R3X]": "[N;H1,H0;+0:2]",
            },
        ],
        "mechanism": [
            {
                "template": "[R1X].[C:10]-[O:3]-[C;H0;+0:4](=[O;+0:5])-[O:6]-[C;H0;+0:7](=[O;+0:8])-[O;H0;+0:9]-[C:11]>>[R2X]-[C;H0;+0:4](=[O;+0:5])-[O;+0:3]-[C:10].[O;-1:6]-[C;H0;+0:7](=[O;+0:8])-[O;H0;+0:9]-[C:11]",
                "description": "Amine reacted with dicarbonate",
            },
            {
                "template": "[O;-1:6]-[C;H0;+0:7](=[O;+0:8])-[O;H0;+0:9]-[C:11]>>[O;+0:6]=[C;H0;+0:7]=[O;+0:8].[O;H0;-1:9]-[C:11]",
                "description": "Carbon dioxide released",
            },
            {
                "template": "[R2X]-[C;H0;+0:4](=[O;+0:5])-[O;+0:3]-[C:10].[O;H0;-1:9]-[C:11]>>[R3X]-[C;H0;+0:4](=[O;+0:5])-[O;+0:3]-[C:10].[O;H1;+0:9]-[C:11]",
                "description": "Proton exchange",
            },
        ],
    },
    "Amide from ester + amine": {
        "scope": [
            {
                "name": "any amine + ester",
                "[R5X]": "[O;H0;+0:1]",
                "[R6X]": "[O;H0;-1:1]",
            },
        ],
        "mechanism": [
            {
                "template": "[#7;H1,H2,H3;+0:2].[R5X]-[C;H0;+0:14](-[#6;+0:13])=[O;H0;+0:15]>>[#6;+0:13]-[C;H0;+0:14](-[O;H0;-1:15])(-[R5X])-[#7;H1,H2,H3;+1:2]",
                "description": "Addition of amine to carbonyl",
            },
            {
                "template": "[#6;+0:13]-[C;H0;+0:14](-[O;H0;-1:15])(-[R5X])-[#7;H1,H2,H3;+1:2]>>[#6;+0:13]-[C;H0;+0:14](-[O;H0;-1:15])(-[R5X])-[#7;H0,H1,H2;+0:2]",
                "description": "Deprotonation of amine",
            },
            {
                "template": "[#6;+0:13]-[C;H0;+0:14](-[O;H0;-1:15])(-[R5X])-[#7;H0,H1,H2;+0:2]>>([#6;+0:13]-[C;H0;+0:14](=[O;H0;+0:15])-[#7;H0,H1,H2;+0:2].[R6X])",
                "description": "Alkoxide leaving",
            },
        ],
    },
    "Isocyanate additions": {
        "scope": [
            {
                "name": "amines",
                "[R1X]": "[#7;H3,H2,H1;+0:4]",
                "[R2X]": "[#7;H3,H2,H1;+1:4]",
                "[R3X]": "[#7;H2,H1,H0;+0:4]",
            },
            {
                "name": "water/alcohol",
                "[R1X]": "[O;H2,H1;+0:4]",
                "[R2X]": "[O;H2,H1;+1:4]",
                "[R3X]": "[O;H1,H0;+0:4]",
            },
        ],
        "mechanism": [
            {
                "template": "[N;H1,H0;+0:1]=[C;H0;+0:2]=[O;H0;+0:3].[R1X]>>[N;H1,H0;-1:1]-[C;H0;+0:2](-[R2X])=[O;H0;+0:3]",
                "description": "Amine addition to isocyanate",
            },
            {
                "template": "[N;H1,H0;-1:1]-[C;H0;+0:2](-[R2X])=[O;H0;+0:3]>>[N;H2,H1;+0:1]-[C;H0;+0:2](-[R3X])=[O;H0;+0:3]",
                "description": "Proton exchange",
            },
        ],
    },
    "Deprotection (N-Cbz) ": {
        "scope": [
            {
                "name": "amines",
                "[R1X]": "[#7;H0,H1;+0:6]",
                "[R2X]": "[#7;H0,H1;-1:6]",
                "[R3X]": "[#7;H1,H2;+0:6]",
            },
        ],
        "mechanism": [
            {
                "template": "[O;H0;+0:1]=[C;H0;+0:2](-[R1X])-[O;H0;+0:3]-[C;H2;+0:5].[Cl;H1;+0:4]>>[O;H1;+1:1]=[C;H0;+0:2](-[R1X])-[O;H0;+0:3]-[C;H2;+0:5].[Cl;H0;-1:4]",
                "description": "Proton exchange",
            },
            {
                "template": "[O;H1;+1:1]=[C;H0;+0:2](-[R1X])-[O;H0;+0:3]-[C;H2;+0:5].[Cl;H0;-1:4]>>([O;H1;+0:1]-[C;H0;+0:2](-[R1X])=[O;H0;+0:3].[C;H2;+0:5]-[Cl;H0;+0:4])",
                "description": "Nucleophilic substitution",
            },
            {
                "template": "[O;H1;+0:1]-[C;H0;+0:2](-[R1X])=[O;H0;+0:3]>>[O;H0;-1:1]-[C;H0;+0:2](-[R1X])=[O;H0;+0:3]",
                "description": "Deprotonation of carboxylic acid",
            },
            {
                "template": "[O;H0;-1:1]-[C;H0;+0:2](-[R1X])=[O;H0;+0:3]>>[R2X].[O;H0;+0:1]=[C;H0;+0:2]=[O;H0;+0:3]",
                "description": "CO2 evolution",
            },
            {"template": "[R2X]>>[R3X]", "description": "Protonation of amide anion"},
        ],
    },
    "Mitsunobu Reaction": {
        "scope": [
            {
                "name": "carboxylic acid and alcohol",
                "[R1X]": "[C;H0;+0:5](=[O;H0;+0:8])",
                "[R2X]": "[C;H2,H1,H0;+0:6]",
            },
            {
                "name": "phenol and alcohol",
                "[R1X]": "[c;H0;+0:5]",
                "[R2X]": "[C;H2,H1,H0;+0:6]",
            },
        ],
        "mechanism": [
            {
                "template": "[P;H0;+0:1].[N;H0;+0:2]=[N;H0;+0:3]>>[P;H0;+1:1]-[N;H0;+0:2]-[N;H0;-1:3]",
                "description": "Reaction between phosphine and diethylazodicarboxylate (DEAD)",
            },
            {
                "template": "[P;H0;+1:1]-[N;H0;+0:2]-[N;H0;-1:3].[R1X]-[O;H1;+0:4]>>[P;H0;+1:1]-[N;H0;+0:2]-[N;H1;+0:3].[R1X]-[O;H0;-1:4]",
                "description": "Proton exchange",
            },
            {
                "template": "[P;H0;+1:1]-[N;H0;+0:2]-[N;H1;+0:3].[R2X]-[O;H1;+0:7]>>[R2X]-[O;H1;+1:7]-[P;H0;+0:1]-[N;H0;+0:2]-[N;H1;+0:3]",
                "description": "Forming a complex with the alcohol",
            },
            {
                "template": "[R2X]-[O;H1;+1:7]-[P;H0;+0:1]-[N;H0;+0:2]-[N;H1;+0:3]>>[R2X]-[O;H0;+0:7]-[P;H0;+0:1]-[N;H1;+1:2]-[N;H1;+0:3]",
                "description": "Proton exchange",
            },
            {
                "template": "[R2X]-[O;H0;+0:7]-[P;H0;+0:1]-[N;H1;+1:2]-[N;H1;+0:3]>>[R2X]-[O;H0;+1:7]=[P;H0;+0:1].[N;H1;+0:2]-[N;H1;+0:3]",
                "description": "Bond cleavage",
            },
            {
                "template": "[R2X]-[O;H0;+1:7]=[P;H0;+0:1].[R1X]-[O;H0;-1:4]>>[O;H0;+0:7]=[P;H0;+0:1].[R1X]-[O;H0;+0:4]-[R2X]",
                "description": "SN2 type reaction",
            },
        ],
    },
    "Friedel-Crafts Acylation": {
        "scope": [
            {"name": "acyl chloride", "[R1X]": "[C;H3,H2,H1,H0;+0:18]"},
        ],
        "mechanism": [
            {
                "template": "[Al;+0:1](-[Cl:2])(-[Cl:3])-[Cl:4].[Cl;H0;+0:5]-[C;H0;+0:6](=[O;H0;+0:7])[R1X]>>[O;H0;+0:7]=[C;H0;+0:6]([R1X])-[Cl;+1:5]-[Al;-1:1](-[Cl:2])(-[Cl:3])-[Cl:4]",
                "description": "Lewis acid activation",
            },
            {
                "template": "[O;H0;+0:7]=[C;H0;+0:6]([R1X])-[Cl;+1:5]-[Al;-1:1](-[Cl:2])(-[Cl:3])-[Cl:4]>>[O;H0;+0:7]=[C;H0;+1:6]([R1X]).[Cl;+0:5]-[Al;-1:1](-[Cl:2])(-[Cl:3])-[Cl:4]",
                "description": "AlCl4 leaves",
            },
            {
                "template": "[O;H0;+0:7]=[C;H0;+1:6]([R1X]).[c;H1;+0:8]:[c;H1;+0:9]>>[c;H1;+1:8]:[c;H1;+0:9]-[C;H0;+0:6](=[O;H0;+0:7])[R1X]",
                "description": "Addition to aromatic ring",
            },
            {
                "template": "[c;H1;+1:8]:[c;H1;+0:9]-[C;H0;+0:6](=[O;H0;+0:7])[R1X].[Cl;+0:5]-[Al;-1:1](-[Cl:2])(-[Cl:3])-[Cl:4]>>[c;H1;+0:8]:[c;H0;+0:9]-[C;H0;+0:6](=[O;H0;+0:7])[R1X].[Cl;H1;+0:5].[Al;+0:1](-[Cl:2])(-[Cl:3])-[Cl:4]",
                "description": "Proton exchange",
            },
        ],
    },
    "Condensation (Aldol)": {
        "scope": [
            {"name": "ketone", "[R1X]": "[C;+0:6]"},
            {"name": "ester", "[R1X]": "[O;H0;+0:6]"},
        ],
        "mechanism": [
            {
                "template": "[O;H0;+0:1]=[C;H0;+0:2](-[R1X])-[C;H2;+0:3]>>[O;H0;+0:1]=[C;H0;+0:2](-[R1X])-[C;H1;-1:3]",
                "description": "Deprotonation of alpha position carbon",
            },
            {
                "template": "[O;H0;+0:1]=[C;H0;+0:2](-[R1X])-[C;H1;-1:3].[C;!$(C(=O)O);+0:4]=[O;H0;+0:5]>>[O;H0;+0:1]=[C;H0;+0:2](-[R1X])-[C;H1;+0:3]-[C;!$(C(=O)O);+0:4]-[O;H0;-1:5]",
                "description": "Addition of enolate",
            },
            {
                "template": "[O;H0;+0:1]=[C;H0;+0:2](-[R1X])-[C;H1;+0:3]-[C;!$(C(=O)O);+0:4]-[O;H0;-1:5]>>[O;H0;+0:1]=[C;H0;+0:2](-[R1X])-[C;H1;+0:3]-[C;!$(C(=O)O);+0:4]-[O;H1;+0:5]",
                "description": "Protonation of alkoxide 1",
            },
            {
                "template": "[O;H0;+0:1]=[C;H0;+0:2](-[R1X])-[C;H1;+0:3]-[C;!$(C(=O)O);+0:4]-[O;H1;+0:5]>>[O;H0;+0:1]=[C;H0;+0:2](-[R1X])-[C;H0;-1:3]-[C;!$(C(=O)O);+0:4]-[O;H2;+1:5]",
                "description": "Protonation of alkoxide 2",
            },
            {
                "template": "[O;H0;+0:1]=[C;H0;+0:2](-[R1X])-[C;H0;-1:3]-[C;!$(C(=O)O);+0:4]-[O;H2;+1:5]>>[O;H0;+0:1]=[C;H0;+0:2](-[R1X])-[C;H0;+0:3]=[C;!$(C(=O)O);+0:4].[O;H2;+0:5]",
                "description": "Condensation",
            },
        ],
    },
    "Condensation (Enamine)": {
        "scope": [
            {"name": "ketone", "[R1X]": "[C;H0;+0:2]"},
            {"name": "aldehyde", "[R1X]": "[C;H1;+0:2]"},
        ],
        "mechanism": [
            {
                "template": "[O;H0;+0:1]=[R1X]-[C;H2;+0:3].[N;H1;+0:4]>>[O;H0;-1:1]-[R1X](-[N;H1;+1:4])-[C;H2;+0:3]",
                "description": "addition of secondary amine into ketone",
            },
            {
                "template": "[O;H0;-1:1]-[R1X](-[N;H1;+1:4])-[C;H2;+0:3]>>[O;H1;+0:1]-[R1X](-[N;H0;+0:4])-[C;H2;+0:3]",
                "description": "proton transfer",
            },
            {
                "template": "[O;H1;+0:1]-[R1X](-[N;H0;+0:4])-[C;H2;+0:3]>>[O;H2;+0:1].[R1X](-[N;H0;+0:4])=[C;H1;+0:3]",
                "description": "condensation",
            },
            {
                "template": "[R1X](-[N;H0;+0:4])=[C;H1;+0:3]>>[R1X](=[N;H0;+1:4])-[C;H1;-1:3]",
                "description": "resonance",
            },
        ],
    },
    "Hydrolysis": {
        "scope": [
            {
                "name": "halide",
                "[R1X]": "[Cl,Br,I;H0;+0:5]",
                "[R2X]": "[Cl,Br,I;H0;-1:5]",
            },
            {"name": "water", "[R1X]": "[O;H2;+1:5]", "[R2X]": "[O;H2;+0:5]"},
            {"name": "ester", "[R1X]": "[O;H0;+0:5]", "[R2X]": "[O;H0;-1:5]"},
            {
                "name": "amides",
                "[R1X]": "[#7;H2,H1,H0;+0:5]",
                "[R2X]": "[#7;H2,H1,H0;-1:5]",
            },
        ],
        "mechanism": [
            {
                "template": "[O;H0;+0:2]=[C;+0:1].[O;H1;-1:3]>>[O;H0;-1:2][C;+0:1][O;H1;+0:3]",
                "description": "hydroxide attack of carbonyl",
            },
            {
                "template": "[O;H0;-1:2]-[C;+0:1]-[R1X]>>([O;H0;+0:2]=[C;+0:1].[R2X])",
                "description": "Leaving or workup",
            }
        ],
    },
    "Wittig Olefination": {
        "scope": [
            {
                "name": "methylide,ethylide",
                "[R1X]": "[C;H2,H1;+0:2]",
                "[R2X]": "[C;H1,H0;-1:2]",
                "[R3X]": "[C;H1,H0;+0:2]",
            }
        ],
        "mechanism": [
            {
                "template": "[P;+1:1]-[R1X][C:5]>>[P;+1:1]-[R2X][C:5]",
                "description": "deprotonation of alpha carbon",
            },
            {
                "template": "[P;+1:1]-[R2X][C:5]>>[P;+0:1]=[R3X][C:5]",
                "description": "tautomerization of ylide",
            },
            {
                "template": "[P;+0:1]=[R3X][C:5].[C;+0:3]=[O;H0;+0:4]>>[P;+0:1]1-[R3X]([C:5])-[C;+0:3]-[O;H0;+0:4]1",
                "description": "Cycloaddition",
            },
            {
                "template": "[P;+0:1]1-[R3X]([C:5])-[C;+0:3]-[O;H0;+0:4]1>>[P;+0:1]=[O;H0;+0:4].[C:5][R3X]=[C;+0:3]",
                "description": "Cyclo-reversion",
            },
        ],
    },
    # "Glorius Saturation": {
    #     "scope": [
    #         {
    #             "name": "benzene",
    #             "[R1X]": "[c;H1;+0:2]",
    #             "[R2X]": "[C;H2;+0:2]",
    #         },
    #         {
    #             "name": "pyridine",
    #             "[R1X]": "[n;+0:2]",
    #             "[R2X]": "[N;H1;+0:2]",
    #         }
    #     ],
    #     "mechanism": [
    #         {
    #             "template": "[c;+0:1]1[R1X][c;H1;+0:3][c;+0:6][c;H1;+0:4][c;H1;+0:5]1>>[C;+0:1]1[R2X][C;H2;+0:3][C;H2;+0:4][C;H2;+0:5][C;+0:6]1",
    #             "description": "saturation",
    #         },
    #     ],
    # },
    "Pinner Reaction": {
        "scope": [{"name": "alkyl alcohol", "[R1X]": "[C;+0:3]"}],
        "mechanism": [
            {
                "template": "[C;H0;+0:1]#[N;H0;+0:2]>>[C;H0;+0:1]#[N;H1;+1:2]",
                "description": "Protonation of nitrile",
            },
            {
                "template": "[C;H0;+0:1]#[N;H1;+1:2].[R1X]-[O;H1;+0:4]>>[R1X]-[O;H1;+1:4]-[C;H0;+0:1]=[N;H1;+0:2]",
                "description": "Addition of alcohol or amine",
            },
            {
                "template": "[R1X]-[O;H1;+1:4]-[C;H0;+0:1]=[N;H1;+0:2]>>[R1X]-[O;H0;+0:4]-[C;H0;+0:1]=[N;H1;+0:2]",
                "description": "Deprotonation or proton exchange",
            },
            {
                "template": "[R1X]-[O;H0;+0:4]-[C;H0;+0:1]=[N;H1;+0:2]>>[R1X]-[O;H0;+0:4]-[C;H0;+0:1]=[N;H2;+1:2]",
                "description": "protonation into Pinner salt",
            },
        ],
    },
    "Pinner Additions": {
        "scope": [{"name": "alkyl alcohol", "[R1X]": "[C;+0:10]"}],
        "mechanism": [
            {
                "template": "[R1X]-[O;H0;+0:4]-[C;+0:1]=[N;H2;+1:2].[N;+0:3]>>[R1X]-[O;H0;-1:4].[N;+1:3][C;+0:1]=[N;H2;+1:2]",
                "description": "amine addition into Pinner salt",
            },
            {
                "template": "[R1X]-[O;H0;+0:4]-[C;+0:1]=[N;H2;+1:2].[O;H2;+0:3]>>[R1X]-[O;H0;+0:4]-[C;+0:1]=[O;H0;+0:3].[N;H3;+0:2]",
                "description": "water addition into Pinner salt",
            },
        ],
    },
    "Pictet-Spengler": {
        "scope": [{"name": "alkyl alcohol", "[R1X]": "[C;+0:10]"}],
        "mechanism": [
            {
                "template": "[C:2](-,:[C:3]-,:[#6;+0:4]1=,:[#6;H1:5]-,:[#7;+0:6]-,:[#6:11]2=,:[#6:12]-,:1-,:[#6:13]=,:[#6:14]-,:[#6:15]=,:[#6:16]-,:2)-[N;H2,H1;+0:1].[C;H1;+0:7]=[O;H0;+0:8]>>[C:2](-[C:3]-[#6;+0:4]1=,:[#6;H1:5]-,:[#7;+0:6]-,:[#6:11]2=,:[#6:12]-,:1-,:[#6:13]=,:[#6:14]-,:[#6:15]=,:[#6:16]-,:2)-[N;H1,H0;+1:1]=[C;H1;+0:7].[O;H2;+0:8]",
                "description": "iminium formation",
            },
            {
                "template": "[C:2](-[C:3]-[#6;+0:4]1=,:[#6;H1:5]-,:[#7;+0:6]-,:[#6:11]2=,:[#6:12]-,:1-,:[#6:13]=,:[#6:14]-,:[#6:15]=,:[#6:16]-,:2)-[N;H1,H0;+1:1]=[C;H1;+0:7]>>[C:2]1-[C:3]-[#6;+0:4]2(-,:[#6;H1:5]=,:[#7&+:6]-,:[#6:11]3=,:[#6:12]-,:2-,:[#6:13]=,:[#6:14]-,:[#6:15]=,:[#6:16]-,:3)-[C&H1&+0:7]-[N;H1,H0;+0:1]-1",
                "description": "spirocycle formation",
            },
            {
                "template": "[C:2]1-[C:3]-[#6;+0:4]2(-,:[#6;H1:5]=,:[#7&+:6]-,:[#6:11]3=,:[#6:12]-,:2-,:[#6:13]=,:[#6:14]-,:[#6:15]=,:[#6:16]-,:3)-[C&H1&+0:7]-[N;H1,H0;+0:1]-1>>[C:2]1-[C:3]-[#6&+:4]2-,:[#6&H1:5](-,:[#7&+0:6]-,:[#6:11]3=,:[#6:12]-,:2-,:[#6:13]=,:[#6:14]-,:[#6:15]=,:[#6:16]-,:3)-[C&H1&+0:7]-[N;H1,H0;+0:1]-1",
                "description": "bond migration",
            },
            {
                "template": "[C:2]1-[C:3]-[#6&+:4]2-,:[#6&H1:5](-,:[#7&+0:6]-,:[#6:11]3=,:[#6:12]-,:2-,:[#6:13]=,:[#6:14]-,:[#6:15]=,:[#6:16]-,:3)-[C&H1&+0:7]-[N;H1,H0;+0:1]-1>>[C:2]1-[C:3]-[#6&+0:4]2=,:[#6&H0:5](-,:[#7&+0:6]-,:[#6:11]3=,:[#6:12]-,:2-,:[#6:13]=,:[#6:14]-,:[#6:15]=,:[#6:16]-,:3)-[C&H1&+0:7]-[N;H1,H0;+0:1]-1",
                "description": "desaturation",
            },
        ],
    },
    "Schotten-Baumann": {
        "scope": [
            {
                "name": "carbonyl/sulfonyl",
                "[R2X]": "[C,S;H0;+0:4]",
            },
        ],
        "mechanism": [
            {
                "template": "[N;H2,H1;+0:1]-[#6;+0:2].[F,Cl,Br,I;H0;+0:3]-[R2X](=[O:5])>>[F,Cl,Br,I;H0;+0:3]-[R2X](-[N;H2,H1;+1:1]-[#6;+0:2])(-[O;-1:5])",
                "description": "Amine addition",
            },
            {
                "template": "[F,Cl,Br,I;H0;+0:3]-[R2X](-[N;H2,H1;+1:1]-[#6;+0:2])(-[O;-1:5])>>[F,Cl,Br,I;H0;+0:3]-[R2X](-[N;H1,H0;+0:1]-[#6;+0:2])(-[O;-1:5])",
                "description": "Deprotonation",
            },
            {
                "template": "[F,Cl,Br,I;H0;+0:3]-[C,S;H0;+0:4](-[N;H1,H0;+0:1]-[#6;+0:2])(-[O;-1:5])>>[F,Cl,Br,I;H0;-1:3].[C,S;H0;+0:4](-[N;H1,H0;+0:1]-[#6;+0:2])(=[O;+0:5])",
                "description": "Halide leaves",
            },
        ],
    },
    "Keto alpha-alkylation": {
        "scope": [
            {
                "name": "alkyl halide",
                "[R2X]": "[C;+0:5]",
            },
        ],
        "mechanism": [
            {
                "template": "[O;H0;+0:1]=[C;H0;+0:2](-[C;H3,H2,H1;+0:3])-[C;+0:6]>>[O;H0;+0:1]=[C;H0;+0:2](-[C;H2,H1,H0;-1:3])-[C;+0:6]",
                "description": "Deprotonation of alpha position carbon",
            },
            {
                "template": "[O;H0;+0:1]=[C;H0;+0:2](-[C;H2,H1,H0;-1:3])-[C;+0:6]>>[O;H0;-1:1]-[C;H0;+0:2](=[C;H2,H1,H0;+0:3])-[C;+0:6]",
                "description": "tautomerization",
            },
            {
                "template": "[O;H0;-1:1]-[C;H0;+0:2](=[C;H2,H1,H0;+0:3])-[C;+0:6].[Cl,Br,I;H0;+0:4]-[R2X]>>[O;H0;+0:1]=[C;H0;+0:2](-[C;H2,H1,H0;+0:3]-[R2X])-[C;+0:6].[Cl,Br,I;H0;-1:4]",
                "description": "Substitution",
            },
        ],
    },
    "Deprotection (ester -> acid)": {
        "scope": [
            {
                "name": "acid alpha CH2 - Me deprotection",
                "[R1X]": "[C;H2;+0:7]",
                "[R2X]": "[C;H3;+0:3]",
            },
            {
                "name": "acid alpha CH3 - Me deprotection",
                "[R1X]": "[C;H3;+0:7]",
                "[R2X]": "[C;H3;+0:3]",
            },
        ],
        "mechanism": [
            {
                "template": "[Li,Na,K:1][O;H1;+0:2].[R2X]-[O;H0;+0:4]-[C;H0;+0:5](-[R1X])=[O;H0;+0:6]>>[Li,Na,K;+1:1].[R2X]-[O;H0;+0:4]-[C;H0;+0:5](-[R1X])(-[O;H1;+0:2])-[O;H0;-1:6]",
                "description": "Hydroxide addition",
            },
            {
                "template": "[R2X]-[O;H0;+0:4]-[C;H0;+0:5](-[R1X])(-[O;H1;+0:2])-[O;H0;-1:6]>>[R2X]-[O;H0;-1:4].[C;H0;+0:5](-[R1X])(-[O;H1;+0:2])=[O;H0;+0:6]",
                "description": "Alkoxide leaves",
            },
        ],
    },
    "Mannich Reaction": {
        "scope": [
            {
                "name": "alkyl/aryl ketone and amine, nonenolizable aldehyde carbonyl",
                "[R1X]": "[#6;+0:4]",
                "[R2X]": "[#6;+0:9]",
            },
        ],
        "mechanism": [
            {
                "template": "[C;H2,H1;!$(C(=O)O);!$(C[C;H2,H1]);+0:2]=[O;H0;+0:3].[R1X]-[N;H2,H1;+0:5]>>[C;H2,H1;!$(C(=O)O);!$(C[C;H2,H1]);+0:2](-[O;H0;-1:3])-[N;H2,H1;+1:5]-[R1X]",
                "description": "Amine addition to carbonyl group",
            },
            {
                "template": "[C;H2,H1;!$(C(=O)O);!$(C[C;H2,H1]);+0:2](-[O;H0;-1:3])-[N;H2,H1;+1:5]-[R1X]>>[C;H2,H1;!$(C(=O)O);!$(C[C;H2,H1]);+0:2](-[O;H1;+0:3])-[N;H1,H0;+0:5]-[R1X]",
                "description": "loss of hydrogen from protonated amine and protonation of alcohol",
            },
            {
                "template": "[C;H2,H1;!$(C(=O)O);!$(C[C;H2,H1]);+0:2](-[O;H1;+0:3])-[N;H1,H0;+0:5]-[R1X]>>[C;H2,H1;!$(C(=O)O);!$(C[C;H2,H1]);+0:2](-[O;H2;+1:3])-[N;H1,H0;+0:5]-[R1X]",
                "description": "protonation of alcohol",
            },
            {
                "template": "[C;H2,H1;!$(C(=O)O);!$(C[C;H2,H1]);+0:2](-[O;H2;+1:3])-[N;H1,H0;+0:5]-[R1X]>>[C;H2,H1;!$(C(=O)O);!$(C[C;H2,H1]);+0:2]=[N;H1,H0;+1:5]-[R1X].[O;H2;+0:3]",
                "description": "loss of water ",
            },
            {
                "template": "[C;H2,H1;!$(C(=O)O);!$(C[C;H2,H1]);+0:2]=[N;H1,H0;+1:5]-[R1X]>>[C;H2,H1;!$(C(=O)O);!$(C[C;H2,H1]);-1:2]-[N;H1,H0;+0:5]-[R1X]",
                "description": "tautomerization of iminium",
            },
            {
                "template": "[C;H2,H1;!$(C(=O)O);!$(C[C;H2,H1]);-1:2]-[N;H1,H0;+0:5]-[R1X]>>[C;H2,H1;!$(C(=O)O);!$(C[C;H2,H1]);+0:2]=[N;H1,H0;+1:5]-[R1X]",
                "description": "tautomerization of iminium reverse",
            },
            {
                "template": "[C;H3,H2,H1;+0:6]-[C;H0;+0:7](-[R2X])(=[O;H0;+0:8])>>[C;H2,H1,H0;+0:6]=[C;H0;+0:7](-[R2X])(-[O;H0;-1:8])",
                "description": "tautomerization of ketone and loss of proton",
            },
            {
                "template": "[C;H2,H1,H0;+0:6]=[C;H0;+0:7](-[R2X])(-[O;H0;-1:8]).[C;H2,H1;!$(C(=O)O);!$(C[C;H2,H1]);+0:2]=[N;H1,H0;+1:5]-[R1X]>>[C;H2,H1,H0;+0:6](-[C;H2,H1;!$(C(=O)O);!$(C[C;H2,H1]);+0:2]-[N;H1,H0;+0:5]-[R1X])-[C;H0;+0:7](-[R2X])(=[O;H0;+0:8])",
                "description": "formation of mannich product",
            },
        ],
    },
    "Ugi reaction": {
        "scope": [
            {
                "name": "amine",
                "[R1X]": "[#6;+0:4]",
            },
        ],
        "mechanism": [
            {
                "template": "[C;H1;+0:2]=[O;H0;+0:3].[R1X]-[N;H2;+0:5]>>[C;H1;+0:2](-[O;H0;-1:3])-[N;H2;+1:5]-[R1X]",
                "description": "Amine addition to carbonyl group",
            },
            {
                "template": "[C;H1;+0:2](-[O;H0;-1:3])-[N;H2;+1:5]-[R1X]>>[C;H1;+0:2](-[O;H1;+0:3])-[N;H1;+0:5]-[R1X]",
                "description": "loss of hydrogen from protonated amine and protonation of alcohol",
            },
            {
                "template": "[C;H1;+0:2](-[O;H1;+0:3])-[N;H1;+0:5]-[R1X]>>[C;H1;+0:2](-[O;H2;+1:3])-[N;H1;+0:5]-[R1X]",
                "description": "protonation of alcohol ",
            },
            {
                "template": "[C;H1;+0:2](-[O;H2;+1:3])-[N;H1;+0:5]-[R1X]>>[C;H1;+0:2]=[N;H1;+1:5]-[R1X].[O;H2;+0:3]",
                "description": "loss of water ",
            },
            {
                "template": "[C;H1;+0:2]=[N;H1;+1:5]-[R1X]>>[C;H1;+1:2]-[N;H1;+0:5]-[R1X]",
                "description": "tautomerization of iminium",
            },
            {
                "template": "[C;H1;+1:2]-[N;H1;+0:5]-[R1X]>>[C;H1;+0:2]=[N;H1;+1:5]-[R1X]",
                "description": "tautomerization of iminium reverse",
            },
            {
                "template": "[C;H1;+0:2]=[N;H1;+1:5]-[R1X].[N;H0;+1:6]#[C;H0;-1:7]>>[C;H1;+0:2](-[C;H0;+0:7]#[N;H0;+1:6])-[N;H1;+0:5]-[R1X]",
                "description": "formation of nitrilium ion",
            },
            {
                "template": "[C;H1;+0:2](-[C;H0;+0:7]#[N;H0;+1:6])-[N;H1;+0:5]-[R1X].[C;H0;+0:8](=[O;H0;+0:9])-[O;H0;-1:10]>>[C;H1;+0:2](-[C;H0;+0:7](-[O;H0;+0:10]-[C;H0;+0:8](=[O;H0;+0:9]))=[N;H0;+0:6])-[N;H1;+0:5]-[R1X]",
                "description": "formation of ugi product",
            },
            {
                "template": "[C;H1;+0:2](-[C;H0;+0:7](-[O;H0;+0:10]-[C;H0;+0:8](=[O;H0;+0:9]))=[N;H0;+0:6])-[N;H1;+0:5]-[R1X]>>[C;H1;+0:2](-[C;H0;+0:7](=[O;H0;+0:10])-[N;H1;+0:6])-[N;H0;+0:5](-[C;H0;+0:8](=[O;H0;+0:9]))-[R1X]",
                "description": "acyl transfer to final ugi product",
            },
        ],
    },
    "Biginelli reaction": {
        "scope": [{"name": "aryl aldehyde", "[R1X]": "[c;H0;+0:1]"}],
        "mechanism": [
            {
                "template": "[R1X]-[C;H1;+0:2]=[O;H0;+0:3].[N;H2;+0:4]-[C;H0+0:5](=[O;H0;+0:6])-[N;H2;+0:7]>>[R1X]-[C;H1;+0:2](-[N;H2;+1:4]-[C;H0;+0:5](=[O;H0;+0:6])-[N;H2;+0:7])-[O;H0;-1:3]",
                "description": "Amine addition to carbonyl group",
            },
            {
                "template": "[R1X]-[C;H1;+0:2](-[N;H2;+1:4]-[C;H0;+0:5](=[O;H0;+0:6])-[N;H2;+0:7])-[O;H0;-1:3]>>[R1X]-[C;H1;+0:2](-[N;H1;+0:4]-[C;H0;+0:5](=[O;H0;+0:6])-[N;H2;+0:7])-[O;H1;+0:3]",
                "description": "loss of hydrogen from protonated amine and protonation of alcohol",
            },
            {
                "template": "[R1X]-[C;H1;+0:2](-[N;H1;+0:4]-[C;H0;+0:5](=[O;H0;+0:6])-[N;H2;+0:7])-[O;H1;+0:3]>>[R1X]-[C;H1;+0:2](-[N;H1;+0:4]-[C;H0;+0:5](=[O;H0;+0:6])-[N;H2;+0:7])-[O;H2;+1:3]",
                "description": "protonation of alcohol ",
            },
            {
                "template": "[R1X]-[C;H1;+0:2](-[N;H1;+0:4]-[C;H0;+0:5](=[O;H0;+0:6])-[N;H2;+0:7])-[O;H2;+1:3]>>[R1X]-[C;H1;+0:2](=[N;H1;+1:4]-[C;H0;+0:5](=[O;H0;+0:6])-[N;H2;+0:7]).[O;H2;+0:3]",
                "description": "loss of water ",
            },
            {
                "template": "[O;H0;+0:8]=[C:9](-[C:10])-[C;H2;+0:11]-[C:12](-[O:13]-[C:14]-[C:15])=[O:16]>>[O;H0;+0:8]=[C:9](-[C:10])-[CH1;-1:11]-[C:12](-[O:13]-[C:14]-[C:15])=[O:16]",
                "description": "formation of ketoester enol",
            },
            {
                "template": "[O;H0;+0:8]=[C:9](-[C:10])-[CH1;-1:11]-[C:12](-[O:13]-[C:14]-[C:15])=[O:16].[R1X]-[C;H1;+0:2](=[N;H1;+1:4]-[C;H0;+0:5](=[O;H0;+0:6])-[N;H2;+0:7])>>[R1X]-[C;H1;+0:2](-[C;H1;+0:11](-[C:9](=[O;H0;+0:8])(-[C:10]))(-[C:12](-[O:13]-[C:14]-[C:15])=[O:16]))(-[N;H1;+0:4]-[C;H0;+0:5](=[O;H0;+0:6])-[N;H2;+0:7])",
                "description": "addition into iminium",
            },
            {
                "template": "[R1X]-[C;H1;+0:2](-[C;H1;+0:11](-[C:9](=[O;H0;+0:8])(-[C:10]))(-[C:12](-[O:13]-[C:14]-[C:15])=[O:16]))(-[N;H1;+0:4]-[C;H0;+0:5](=[O;H0;+0:6])-[N;H2;+0:7])>>[R1X]-[C;H1;+0:2](-[C;H1;+0:11](-[C:9]1(-[O;H0;-1:8])(-[C:10]))(-[C:12](-[O:13]-[C:14]-[C:15])=[O:16]))(-[N;H1;+0:4]-[C;H0;+0:5](=[O;H0;+0:6])-[N;H2;+1:7]1)",
                "description": "amine addition into ketone cyclization",
            },
            {
                "template": "[R1X]-[C;H1;+0:2](-[C;H1;+0:11](-[C:9]1(-[O;H0;-1:8])(-[C:10]))(-[C:12](-[O:13]-[C:14]-[C:15])=[O:16]))(-[N;H1;+0:4]-[C;H0;+0:5](=[O;H0;+0:6])-[N;H2;+1:7]1)>>[R1X]-[C;H1;+0:2](-[C;H1;+0:11](-[C:9]1(-[O;H1;+0:8])(-[C:10]))(-[C:12](-[O:13]-[C:14]-[C:15])=[O:16]))(-[N;H1;+0:4]-[C;H0;+0:5](=[O;H0;+0:6])-[N;H1;+0:7]1)",
                "description": "protonation of alcohol and loss of amine proton",
            },
            {
                "template": "[R1X]-[C;H1;+0:2](-[C;H1;+0:11](-[C:9]1(-[O;H1;+0:8])(-[C:10]))(-[C:12](-[O:13]-[C:14]-[C:15])=[O:16]))(-[N;H1;+0:4]-[C;H0;+0:5](=[O;H0;+0:6])-[N;H1;+0:7]1)>>[R1X]-[C;H1;+0:2](-[C;H1;+0:11](-[C:9]1(-[O;H2;+1:8])(-[C:10]))(-[C:12](-[O:13]-[C:14]-[C:15])=[O:16]))(-[N;H1;+0:4]-[C;H0;+0:5](=[O;H0;+0:6])-[N;H1;+0:7]1)",
                "description": "protonation of alcohol",
            },
            {
                "template": "[R1X]-[C;H1;+0:2](-[C;H1;+0:11](-[C:9]1(-[O;H2;+1:8])(-[C:10]))(-[C:12](-[O:13]-[C:14]-[C:15])=[O:16]))(-[N;H1;+0:4]-[C;H0;+0:5](=[O;H0;+0:6])-[N;H1;+0:7]1)>>[R1X]-[C;H1;+0:2](-[C;H0;+0:11](=[C:9]1(-[C:10]))(-[C:12](-[O:13]-[C:14]-[C:15])=[O:16]))(-[N;H1;+0:4]-[C;H0;+0:5](=[O;H0;+0:6])-[N;H1;+0:7]1).[O;H2;+0:8]",
                "description": "loss of water and formation of Biginelli product",
            },
        ],
    },
    "Grieco Coupling": {
        "scope": [{"name": "any aldehyde", "[R1X]": "[#6;+0:1]"}],
        "mechanism": [
            {
                "template": "[R1X]-[C;H1;+0:2]=[O;H0;+0:3].[N;H2+0:4]-[c;H0;+0:5]:[c;H1+0:6]>>[R1X]-[C;H1;+0:2](-[N;H2;+1:4]-[c;H0;+0:5]:[c;H1;+0:6])-[O;H0;-1:3]",
                "description": "Amine addition to carbonyl group",
            },
            {
                "template": "[R1X]-[C;H1;+0:2](-[N;H2;+1:4]-[c;H0;+0:5]:[c;H1;+0:6])-[O;H0;-1:3]>>[R1X]-[C;H1;+0:2](-[N;H1;+0:4]-[c;H0;+0:5]:[c;H1;+0:6])-[O;H1;+0:3]",
                "description": "loss of hydrogen from protonated amine and protonation of alcohol",
            },
            {
                "template": "[R1X]-[C;H1;+0:2](-[N;H1;+0:4]-[c;H0;+0:5]:[c;H1;+0:6])-[O;H1;+0:3]>>[R1X]-[C;H1;+0:2](-[N;H1;+0:4]-[c;H0;+0:5]:[c;H1;+0:6])-[O;H2;+1:3]",
                "description": "protonation of alcohol ",
            },
            {
                "template": "[R1X]-[C;H1;+0:2](-[N;H1;+0:4]-[c;H0;+0:5]:[c;H1;+0:6])-[O;H2;+1:3]>>[R1X]-[C;H1;+0:2](=[N;H1;+1:4]-[c;H0;+0:5]:[c;H1;+0:6]).[O;H2;+0:3]",
                "description": "loss of water and formation of immonium ion",
            },
            {
                "template": "[R1X]-[C;H1;+0:2](=[N;H1;+1:4]-[c;H0;+0:5]:[c;H1;+0:6]).[C;+0:7]=[C;+0:8]>>[R1X]-[C;H1;+0:2]1(-[N;H1;+0:4]-[c;H0;+0:5]:[c;H0;+0:6]([C;+0:7]-[C;+0:8]1))",
                "description": "aza-diel alder reaction",
            },
        ],
    },
    "Hantzsch Pyridine Synthesis": {
        "scope": [{"name": "alkyl primary aldehyde", "[R1X]": "[C;H2;+0:17]"}],
        "mechanism": [
            {
                "template": "[C;H0;+0:4](=[O;H0;+0:5])-[C;H2;+0:6]-[C;H0;+0:7](=[O;H0;+0:8])[O;H0;+0:9]>>[C;H0;+0:4](-[O;H0;-1:5])=[C;H1;+0:6][C;H0;+0:7](=[O;H0;+0:8])-[O;H0;+0:9]",
                "description": "deprotonation of betaketoester",
            },
            {
                "template": "[C;H1;+0:2](-[R1X])=[O;H0;+0:3].[C;H0;+0:4](-[O;H0;-1:5])=[C;H1;+0:6][C;H0;+0:7](=[O;H0;+0:8])-[O;H0;+0:9]>>[C;H0;+0:4](=[O;H0;+0:5])-[C;H1;+0:6](-[C;H1;+0:2](-[R1X])-[O;H0;-1:3])-[C;H0;+0:7](=[O;H0;+0:8])-[O;H0;+0:9]",
                "description": "addition into aldehyde",
            },
            {
                "template": "[C;H0;+0:4](=[O;H0;+0:5])-[C;H1;+0:6](-[C;H1;+0:2](-[R1X])-[O;H0;-1:3])-[C;H0;+0:7](=[O;H0;+0:8])-[O;H0;+0:9]>>[C;H0;+0:4](=[O;H0;+0:5])-[C;H1;+0:6](-[C;H1;+0:2](-[R1X])-[O;H1;+0:3])-[C;H0;+0:7](=[O;H0;+0:8])-[O;H0;+0:9]",
                "description": "protonation of alcohol ",
            },
            {
                "template": "[C;H0;+0:4](=[O;H0;+0:5])-[C;H1;+0:6](-[C;H1;+0:2](-[R1X])-[O;H1;+0:3])-[C;H0;+0:7](=[O;H0;+0:8])-[O;H0;+0:9]>>[C;H0;+0:4](=[O;H0;+0:5])-[C;H1;+0:6](-[C;H1;+0:2](-[R1X])-[O;H2;+1:3])-[C;H0;+0:7](=[O;H0;+0:8])-[O;H0;+0:9]",
                "description": "protonation of alcohol",
            },
            {
                "template": "[C;H0;+0:4](=[O;H0;+0:5])-[C;H1;+0:6](-[C;H1;+0:2](-[R1X])-[O;H2;+1:3])-[C;H0;+0:7](=[O;H0;+0:8])-[O;H0;+0:9]>>[C;H0;+0:4](=[O;H0;+0:5])-[C;H0;+0:6](=[C;H1;+0:2](-[R1X]))-[C;H0;+0:7](=[O;H0;+0:8])-[O;H0;+0:9].[O;H2;+0:3]",
                "description": "loss of water and formation of condensation product",
            },
            {
                "template": "[C;H0;+0:11](=[O;H0;+0:12])-[C;H2;+0:13]-[C;H0;+0:14](=[O;H0;+0:15])-[O;H0;+0:16].[N;H3;+0:10]>>[C;H0;+0:11](-[N;H2;+0:10])(-[O;H1;+0:12])-[C;H2;+0:13]-[C;H0;+0:14](=[O;H0;+0:15])-[O;H0;+0:16]",
                "description": "addition of ammonia into betaketoester",
            },
            {
                "template": "[C;H0;+0:11](-[N;H2;+0:10])(-[O;H1;+0:12])-[C;H2;+0:13]-[C;H0;+0:14](=[O;H0;+0:15])-[O;H0;+0:16]>>[C;H0;+0:11](-[N;H2;+0:10])(-[O;H2;+1:12])-[C;H2;+0:13]-[C;H0;+0:14](=[O;H0;+0:15])-[O;H0;+0:16]",
                "description": "protonation of aminated betaketoester",
            },
            {
                "template": "[C;H0;+0:11](-[N;H2;+0:10])(-[O;H2;+1:12])-[C;H2;+0:13]-[C;H0;+0:14](=[O;H0;+0:15])-[O;H0;+0:16]>>[C;H0;+0:11](-[N;H2;+0:10])=[C;H1;+0:13]-[C;H0;+0:14](=[O;H0;+0:15])-[O;H0;+0:16].[O;H2;+0:12]",
                "description": "loss of water",
            },
            {
                "template": "[C;H0;+0:11](-[N;H2;+0:10])=[C;H1;+0:13]-[C;H0;+0:14](=[O;H0;+0:15])-[O;H0;+0:16].[C;H0;+0:4](=[O;H0;+0:5])-[C;H0;+0:6](=[C;H1;+0:2](-[R1X]))-[C;H0;+0:7](=[O;H0;+0:8])-[O;H0;+0:9]>>[C;H0;+0:11](=[N;H2;+1:10])-[C;H1;+0:13](-[C;H1;+0:2](-[R1X])-[C;H0;+0:6](-[C;H0;+0:7](=[O;H0;+0:8])-[O;H0;+0:9])=[C;H0;+0:4](-[O;H0;-1:5]))-[C;H0;+0:14](=[O;H0;+0:15])-[O;H0;+0:16]",
                "description": "formation of dihydropyridine derivative",
            },
            {
                "template": "[C;H0;+0:11](=[N;H2;+1:10])-[C;H1;+0:13](-[C;H1;+0:2](-[R1X])-[C;H0;+0:6](-[C;H0;+0:7](=[O;H0;+0:8])-[O;H0;+0:9])=[C;H0;+0:4](-[O;H0;-1:5]))-[C;H0;+0:14](=[O;H0;+0:15])-[O;H0;+0:16]>>[C;H0;+0:11](=[N;H1;+0:10])-[C;H1;+0:13](-[C;H1;+0:2](-[R1X])-[C;H0;+0:6](-[C;H0;+0:7](=[O;H0;+0:8])-[O;H0;+0:9])=[C;H0;+0:4](-[O;H0;-1:5]))-[C;H0;+0:14](=[O;H0;+0:15])-[O;H0;+0:16]",
                "description": "deprotonation of amine",
            },
            {
                "template": "[C;H0;+0:11](=[N;H1;+0:10])-[C;H1;+0:13](-[C;H1;+0:2](-[R1X])-[C;H0;+0:6](-[C;H0;+0:7](=[O;H0;+0:8])-[O;H0;+0:9])=[C;H0;+0:4](-[O;H0;-1:5]))-[C;H0;+0:14](=[O;H0;+0:15])-[O;H0;+0:16]>>[C;H0;+0:11](=[N;H1;+0:10])-[C;H1;+0:13](-[C;H1;+0:2](-[R1X])-[C;H0;-1:6](-[C;H0;+0:7](=[O;H0;+0:8])-[O;H0;+0:9])-[C;H0;+0:4](=[O;H0;+0:5]))-[C;H0;+0:14](=[O;H0;+0:15])-[O;H0;+0:16]",
                "description": "tautomerization of dihydropyridine derivative",
            },
            {
                "template": "[C;H0;+0:11](=[N;H1;+0:10])-[C;H1;+0:13](-[C;H1;+0:2](-[R1X])-[C;H0;-1:6](-[C;H0;+0:7](=[O;H0;+0:8])-[O;H0;+0:9])-[C;H0;+0:4](=[O;H0;+0:5]))-[C;H0;+0:14](=[O;H0;+0:15])-[O;H0;+0:16]>>[C;H0;+0:11](=[N;H1;+0:10])-[C;H0;-1:13](-[C;H1;+0:2](-[R1X])-[C;H1;+0:6](-[C;H0;+0:7](=[O;H0;+0:8])-[O;H0;+0:9])-[C;H0;+0:4](=[O;H0;+0:5]))-[C;H0;+0:14](=[O;H0;+0:15])-[O;H0;+0:16]",
                "description": "tautomerization of dihydropyridine derivative",
            },
            {
                "template": "[C;H0;+0:11](=[N;H1;+0:10])-[C;H0;-1:13](-[C;H1;+0:2](-[R1X])-[C;H1;+0:6](-[C;H0;+0:7](=[O;H0;+0:8])-[O;H0;+0:9])-[C;H0;+0:4](=[O;H0;+0:5]))-[C;H0;+0:14](=[O;H0;+0:15])-[O;H0;+0:16]>>[C;H0;+0:11](-[N;H1;+0:10]1)=[C;H0;+0:13](-[C;H1;+0:2](-[R1X])-[C;H1;+0:6](-[C;H0;+0:7](=[O;H0;+0:8])-[O;H0;+0:9])-[C;H0;+0:4]1(-[O;H1;+0:5]))-[C;H0;+0:14](=[O;H0;+0:15])-[O;H0;+0:16]",
                "description": "cyclization",
            },
            {
                "template": "[C;H0;+0:11](-[N;H1;+0:10]1)=[C;H0;+0:13](-[C;H1;+0:2](-[R1X])-[C;H1;+0:6](-[C;H0;+0:7](=[O;H0;+0:8])-[O;H0;+0:9])-[C;H0;+0:4]1(-[O;H1;+0:5]))-[C;H0;+0:14](=[O;H0;+0:15])-[O;H0;+0:16]>>[C;H0;+0:11](-[N;H1;+0:10]1)=[C;H0;+0:13](-[C;H1;+0:2](-[R1X])-[C;H0;+0:6](-[C;H0;+0:7](=[O;H0;+0:8])-[O;H0;+0:9])=[C;H0;+0:4]1)-[C;H0;+0:14](=[O;H0;+0:15])-[O;H0;+0:16].[O;H2;+0:5]",
                "description": "loss of water and formation of pyridine product",
            },
        ],
    },
    "Kabachnik-Fields-anilines": {
        "scope": [
            {
                "name": "aniline - ketone",
                "[R1X]": "[c;H0;+0:9]",
                "[R2X]": "[C;H0;+0;!$(C(=O)O):2]",
            },
            {
                "name": "aniline - aldehyde",
                "[R1X]": "[c;H0;+0:9]",
                "[R2X]": "[C;H1;+0:2]",
            },
        ],
        "mechanism": [
            {
                "template": "[R2X]=[O;H0;+0:3].[R1X]-[N;H2,H1;+0:4]>>[R2X](-[N;H2,H1;+1:4]-[R1X])-[O;H0;-1:3]",
                "description": "amine addition to carbonyl group",
            },
            {
                "template": "[R2X](-[N;H2,H1;+1:4]-[R1X])-[O;H0;-1:3]>>[R2X](-[N;H1,H0;+0:4]-[R1X])-[O;H1;+0:3]",
                "description": "protonation of alcohol and deprotonation of amine",
            },
            {
                "template": "[R2X](-[N;H1,H0;+0:4]-[R1X])-[O;H1;+0:3]>>[R2X](-[N;H1,H0;+0:4]-[R1X])-[O;H2;+1:3]",
                "description": "protonation of alcohol",
            },
            {
                "template": "[R2X](-[N;H1,H0;+0:4]-[R1X])-[O;H2;+1:3]>>[R2X](=[N;H1,H0;+1:4]-[R1X]).[O;H2;+0:3]",
                "description": "loss of water and formation of imine",
            },
            {
                "template": "[R2X](=[N;H1,H0;+1:4]-[R1X]).[O;H0;+0:5]-[P;H1;+0:6](=[O;H0;+0:7])-[O;H0;+0:8]>>[O;H0;+0:5]-[P;H0;+0:6](-[R2X](-[N;H1,H0;+0:4]-[R1X]))(=[O;H0;+0:7])-[O;H0;+0:8]",
                "description": "addition of alpha-amino phosphonate to form product",
            },
        ],
    },
    "Kabachnik-Fields-alkylamines": {
        "scope": [
            {"name": "alkyl amines", "[R1X]": "[C;H2,H1,H0;+0:9]"},
        ],
        "mechanism": [
            {
                "template": "[C;H0;+0:2]=[O;H0;+0:3].[O;H0;+0:5]-[P;H1;+0:6](=[O;H0;+0:7])-[O;H0;+0:8]>>[O;H0;+0:5]-[P;H0;+0:6](-[C;H0;+0:2]-[O;H1;+0:3])(=[O;H0;+0:7])-[O;H0;+0:8]",
                "description": "aldehyde addition to alpha-amine phosphonate",
            },
            {
                "template": "[O;H0;+0:5]-[P;H0;+0:6](-[C;H0;+0:2]-[O;H1;+0:3])(=[O;H0;+0:7])-[O;H0;+0:8]>>[O;H0;+0:5]-[P;H0;+0:6](-[C;H0;+0:2]-[O;H2;+1:3])(=[O;H0;+0:7])-[O;H0;+0:8]",
                "description": "protonation of alcohol",
            },
            {
                "template": "[O;H0;+0:5]-[P;H0;+0:6](-[C;H0;+0:2]-[O;H2;+1:3])(=[O;H0;+0:7])-[O;H0;+0:8].[R1X]-[N;H2,H1;+0:4]>>[O;H0;+0:5]-[P;H0;+0:6](-[C;H0;+0:2](-[N;H2,H1;+1:4]-[R1X]))(=[O;H0;+0:7])-[O;H0;+0:8].[O;H2;+0:3]",
                "description": "amine addition and loss of water",
            },
            {
                "template": "[O;H0;+0:5]-[P;H0;+0:6](-[C;H0;+0:2](-[N;H2,H1;+1:4]-[R1X]))(=[O;H0;+0:7])-[O;H0;+0:8]>>[O;H0;+0:5]-[P;H0;+0:6](-[C;H0;+0:2](-[N;H1,H0;+0:4]-[R1X]))(=[O;H0;+0:7])-[O;H0;+0:8]",
                "description": "deprotonation of amine to form product",
            },
        ],
    },
    "Passerini Reaction": {
        "scope": [
            {
                "name": "isocyanide + aldehyde/ketone",
                "[R1X]": "[C;+0:10]",
                "[R2X]": "[C;H1,H0;!$(C(=O)O);+0:2]",
            },
        ],
        "mechanism": [
            {
                "template": "[R2X]=[O;H0;+0:3].[C;H0;+0:4](=[O;H0;+0:5])(-[O;H1;+0:6])>>[R2X]=[O;H0;+1:3]-[O;H0;+0:6]-[C;H0;+0:4](=[O;H0;+0:5])",
                "description": "hydrogen bond formation",
            },
            {
                "template": "[R2X]=[O;H0;+1:3]-[O;H0;+0:6]-[C;H0;+0:4](=[O;H0;+0:5]).[R1X]-[N;H0;+1:8]#[C;H0;-1:9]>>[R2X](-[C;H0;+0:9](-[O;H0;+0:5]-[C;H0;+0:4]=[O;H0;+0:6])=[N;H0;+0:8](-[R1X]))-[O;H1;+0:3]",
                "description": "bond rearrangement",
            },
            {
                "template": "[R2X](-[C;H0;+0:9](-[O;H0;+0:5]-[C;H0;+0:4]=[O;H0;+0:6])=[N;H0;+0:8](-[R1X]))-[O;H1;+0:3]>>[R2X](-[O;H0;+0:3]-[C;H0;+0:4]=[O;H0;+0:6])(-[C;H0;+0:9](=[O;H0;+0:5])-[N;H1;+0:8](-[R1X]))",
                "description": "bond rearrangement to form passerini product ",
            },
        ],
    },
    "Petasis Reaction": {
        "scope": [
            {
                "name": "ArB(OH)2 + amine + aldehyde/ketone",
                "[R1X]": "[c;H0;+0:8]",
                "[R2X]": "[N;H3,H2,H1;+0:4]",
                "[R3X]": "[N;H2,H1,H0;+0:4]",
                "[R4X]": "[N;H2,H1,H0;+1:4]",
                "[R5X]": "[C;H1,H0;!$(C(=O)O);+0:2]",
            },
            {
                "name": "Vinyl-B(OH)2 + amine + aldehyde/ketone",
                "[R1X]": "[C;H1,H0;+0:8]=[C;H2,H1,H0;+0:9]",
                "[R2X]": "[N;H3,H2,H1;+0:4]",
                "[R3X]": "[N;H2,H1,H0;+0:4]",
                "[R4X]": "[N;H2,H1,H0;+1:4]",
                "[R5X]": "[C;H1,H0;!$(C(=O)O);+0:2]",
            },
        ],
        "mechanism": [
            {
                "template": "[R5X]=[O;H0;+0:3].[R2X]>>[R5X](-[O;H1;+0:3])-[R3X]",
                "description": "amine addition into carbonyl group",
            },
            {
                "template": "[R5X](-[O;H1;+0:3])-[R3X].[B;H0;+0:5](-[O:6])(-[O:7])-[R1X]>>[R5X](-[O;H1;+1:3]-[B;H0;-1:5](-[O:6])(-[O:7])-[R1X])-[R3X]",
                "description": "creation of boronate ester",
            },
            {
                "template": "[R5X](-[O;H1;+1:3]-[B;H0;-1:5](-[O:6])(-[O:7])-[R1X])-[R3X]>>([R5X]=[R4X].[O;H1;+0:3]-[B;H0;-1:5](-[O:6])(-[O:7])-[R1X])",
                "description": "cleavage into iminium",
            },
            {
                "template": "[R5X]=[R4X].[O;H1;+0:3]-[B;H0;-1:5](-[O:6])(-[O:7])-[R1X]>>[R5X](-[R1X])-[R3X].[O;H1;+0:3]-[B;H0;+0:5](-[O:6])(-[O:7])",
                "description": "bond abstraction",
            },
        ],
    },
    "Strecker Synthesis": {
        "scope": [
            {
                "name": "CH3 aldehyde (valine) + ammonia",
                "[R1X]": "[C;H3,H2,H1;+0:8]",
                "[R2X]": "[N;H3,H2;+0:4]",
                "[R3X]": "[N;H3,H2;+1:4]",
                "[R4X]": "[N;H2,H1;+0:4]",
                "[R5X]": "[N;H2,H1;+1:4]",
            },
        ],
        "mechanism": [
            {
                "template": "[R1X]-[C;H1;+0:2]=[O;H0;+0:3]>>[R1X]-[C;H1;+0:2]=[O;H1;+1:3]",
                "description": "protonation of aldehyde",
            },
            {
                "template": "[R1X]-[C;H1;+0:2]=[O;H1;+1:3].[R2X]>>[R3X]-[C;H1;+0:2](-[R1X])-[O;H1;+0:3]",
                "description": "addition of amine",
            },
            {
                "template": "[R3X]-[C;H1;+0:2](-[R1X])-[O;H1;+0:3]>>[R4X]-[C;H1;+0:2](-[R1X])-[O;H1;+0:3]",
                "description": "deprotonation of amine ",
            },
            {
                "template": "[R4X]-[C;H1;+0:2](-[R1X])-[O;H1;+0:3]>>[R4X]-[C;H1;+0:2](-[R1X])-[O;H2;+1:3]",
                "description": "protonation of alcohol ",
            },
            {
                "template": "[R4X]-[C;H1;+0:2](-[R1X])-[O;H2;+1:3]>>[R5X]=[C;H1;+0:2](-[R1X]).[O;H2;+0:3]",
                "description": "elimination of water and formation of iminium",
            },
            {
                "template": "[R5X]=[C;H1;+0:2](-[R1X]).[N;H0;+0:5]#[C;H0;-1:6]>>[R4X]-[C;H1;+0:2](-[R1X])-[C;H0;-0:6]#[N;H0;+0:5]",
                "description": "addition of CN",
            },
            {
                "template": "[R4X]-[C;H1;+0:2](-[R1X])-[C;H0;-0:6]#[N;H0;+0:5]>>[R4X]-[C;H1;+0:2](-[R1X])-[C;H0;-0:6]#[N;H1;+1:5]",
                "description": "protonation of nitrile",
            },
            {
                "template": "[R4X]-[C;H1;+0:2](-[R1X])-[C;H0;-0:6]#[N;H1;+1:5].[O;H2;+0:3]>>[R4X]-[C;H1;+0:2](-[R1X])-[C;H0;-0:6](-[O;H2;+1:3])=[N;H1;+0:5]",
                "description": "addition of water",
            },
            {
                "template": "[R4X]-[C;H1;+0:2](-[R1X])-[C;H0;-0:6](-[O;H2;+1:3])=[N;H1;+0:5]>>[R4X]-[C;H1;+0:2](-[R1X])-[C;H0;-0:6](-[O;H1;+0:3])=[N;H1;+0:5]",
                "description": "deprotonation",
            },
            {
                "template": "[R4X]-[C;H1;+0:2](-[R1X])-[C;H0;-0:6](-[O;H1;+0:3])=[N;H1;+0:5]>>[R4X]-[C;H1;+0:2](-[R1X])-[C;H0;-0:6](-[O;H1;+0:3])=[N;H2;+1:5]",
                "description": "protonation",
            },
            {
                "template": "[R4X]-[C;H1;+0:2](-[R1X])-[C;H0;-0:6](-[O;H1;+0:3])=[N;H2;+1:5].[O;H2;+0:7]>>[R4X]-[C;H1;+0:2](-[R1X])-[C;H0;-0:6](-[O;H2;+1:7])(-[O;H1;+0:3])-[N;H2;+0:5]",
                "description": "addition of water",
            },
            {
                "template": "[R4X]-[C;H1;+0:2](-[R1X])-[C;H0;-0:6](-[O;H2;+1:7])(-[O;H1;+0:3])-[N;H2;+0:5]>>[R4X]-[C;H1;+0:2](-[R1X])-[C;H0;-0:6](-[O;H1;+0:7])(-[O;H1;+0:3])-[N;H2;+0:5]",
                "description": "deprotonation",
            },
            {
                "template": "[R4X]-[C;H1;+0:2](-[R1X])-[C;H0;-0:6](-[O;H1;+0:7])(-[O;H1;+0:3])-[N;H2;+0:5]>>[R4X]-[C;H1;+0:2](-[R1X])-[C;H0;-0:6](-[O;H1;+0:7])(-[O;H1;+0:3])-[N;H3;+1:5]",
                "description": "protonation",
            },
            {
                "template": "[R4X]-[C;H1;+0:2](-[R1X])-[C;H0;-0:6](-[O;H1;+0:7])(-[O;H1;+0:3])-[N;H3;+1:5]>>[R4X]-[C;H1;+0:2](-[R1X])-[C;H0;-0:6](=[O;H1;+1:7])(-[O;H1;+0:3]).[N;H3;+0:5]",
                "description": "elimination of amine",
            },
            {
                "template": "[R4X]-[C;H1;+0:2](-[R1X])-[C;H0;-0:6](=[O;H1;+1:7])(-[O;H1;+0:3])>>[R4X]-[C;H1;+0:2](-[R1X])-[C;H0;-0:6](=[O;H0;+0:7])([O;H1;+0:3])",
                "description": "deprotonation to strieker product",
            },
        ],
    },
    "Aldehyde-Alkyne-Amine coupling": {
        "scope": [
            {
                "name": "aldehyde",
                "[R1X]": "[#6;+0:8]",
            },
        ],
        "mechanism": [
            {
                "template": "[C;H1;+0:2](-[R1X])=[O;H0;+0:3]>>[C;H1;+0:2](-[R1X])=[O;H1;+1:3]",
                "description": "protonation of aldehyde",
            },
            {
                "template": "[C;H1;+0:2](-[R1X])=[O;H1;+1:3].[N;H3,H2,H1;+0:4]>>[N;H3,H2,H1;+1:4]-[C;H1;+0:2](-[R1X])-[O;H1;+0:3]",
                "description": "addition of amine",
            },
            {
                "template": "[N;H3,H2,H1;+1:4]-[C;H1;+0:2](-[R1X])-[O;H1;+0:3]>>[N;H2,H1,H0;+0:4]-[C;H1;+0:2](-[R1X])-[O;H1;+0:3]",
                "description": "deprotonation of amine ",
            },
            {
                "template": "[N;H2,H1,H0;+0:4]-[C;H1;+0:2](-[R1X])-[O;H1;+0:3]>>[N;H2,H1,H0;+0:4]-[C;H1;+0:2](-[R1X])-[O;H2;+1:3]",
                "description": "protonation of alcohol ",
            },
            {
                "template": "[N;H2,H1,H0;+0:4]-[C;H1;+0:2](-[R1X])-[O;H2;+1:3]>>[N;H2,H1,H0;+1:4]=[C;H1;+0:2](-[R1X]).[O;H2;+0:3]",
                "description": "elimination of water and formation of iminium",
            },
            {
                "template": "[C;H0;+0:6]#[C;H1;+0:5]>>[C;H0;+0:6]#[C;H0;-1:5]",
                "description": "deprotonation of alkyne",
            },
            {
                "template": "[C;H0;+0:6]#[C;H0;-1:5].[N;H2,H1,H0;+1:4]=[C;H1;+0:2](-[R1X])>>[C;H0;+0:6]#[C;H0;+0:5]-[C;H1:2](-[R1X])-[N;H2,H1,H0;+0:4]",
                "description": "addition of acetylide into iminium",
            },
        ],
    },
    "Favorskii Rearrangement": {
        "scope": [
            {
                "name": "halide",
                "[R1X]": "[Cl,Br,I;H0;+0:4]",
                "[R2X]": "[Cl,Br,I;H0;-1:4]",
            }
        ],
        "mechanism": [
            {
                "template": "[O;H0;+0:1]=[C;H0;+0:2]1-[C;H1;+0:3](-[R1X])-[C;H2;+0:5]-[C;H2;+0:6]-[C;H2;+0:7]-[C;H2;+0:8]1>>[O;H0;-1:1]-[C;H0;+0:2]1-[C;H1;+0:3](-[R1X])-[C;H2;+0:5]-[C;H2;+0:6]-[C;H2;+0:7]-[C;H1;+0:8]=1",
                "description": "beta keto deprotonation",
            },
            {
                "template": "[O;H0;-1:1]-[C;H0;+0:2]1-[C;H1;+0:3](-[R1X])-[C;H2;+0:5]-[C;H2;+0:6]-[C;H2;+0:7]-[C;H1;+0:8]=1>>[O;H0;+0:1]=[C;H0;+0:2]1-[C;H1;+0:3]2-[C;H2;+0:5]-[C;H2;+0:6]-[C;H2;+0:7]-[C;H1;+0:8]-2-1.[R2X]",
                "description": "bridge formation and loss of halide ion",
            },
            {
                "template": "[O;H0;+0:1]=[C;H0;+0:2]1-[C;H1;+0:3]2-[C;H2;+0:5]-[C;H2;+0:6]-[C;H2;+0:7]-[C;H1;+0:8]-2-1.[O;H1;-1:9]>>[O;H0;-1:1]-[C;H0;+0:2](-[O;H1;+0:9])1-[C;H1;+0:3]2-[C;H2;+0:5]-[C;H2;+0:6]-[C;H2;+0:7]-[C;H1;+0:8]-2-1",
                "description": "hydroxylation of ketone",
            },
            {
                "template": "[O;H0;-1:1]-[C;H0;+0:2](-[O;H1;+0:9])1-[C;H1;+0:3]2-[C;H2;+0:5]-[C;H2;+0:6]-[C;H2;+0:7]-[C;H1;+0:8]-2-1>>[O;H0;+0:1]=[C;H0;+0:2](-[O;H1;+0:9])-[C;H1;+0:3]2-[C;H2;+0:5]-[C;H2;+0:6]-[C;H2;+0:7]-[C;H1;-1:8]-2",
                "description": "cyclopropane ring opening",
            },
            {
                "template": "[O;H0;+0:1]=[C;H0;+0:2](-[O;H1;+0:9])-[C;H1;+0:3]2-[C;H2;+0:5]-[C;H2;+0:6]-[C;H2;+0:7]-[C;H1;-1:8]-2>>[O;H0;+0:1]=[C;H0;+0:2](-[O;H1;+0:9])-[C;H1;+0:3]2-[C;H2;+0:5]-[C;H2;+0:6]-[C;H2;+0:7]-[C;H2;+0:8]-2",
                "description": "protonation",
            },
        ],
    },
    "Favorskii Reaction": {
        "scope": [
            {
                "name": "hydroxyl/alkoxide attack",
                "[R1X]": "[O;H1,H0;-1:7]",
                "[R2X]": "[O;H1,H0;+0:7]",
            }
        ],
        "mechanism": [
            {
                "template": "[C;H3,H2,H1;+0:2]-[C&H0&+0:3](-[C;H2,H1;+0:4]-[Cl,Br,I;H0;+0:6])=[O&H0&+0:5]>>[C;H3,H2,H1;+0:2]-[C&H0&+0:3](=[C;H1,H0;+0:4]-[Cl,Br,I;H0;+0:6])-[O&H0&-:5]",
                "description": "alpha-halo ketone tautomerization",
            },
            {
                "template": "[C;H3,H2,H1;+0:2]-[C&H0&+0:3](=[C;H1,H0;+0:4]-[Cl,Br,I;H0;+0:6])-[O&H0&-:5]>>[C;H2,H1,H0;+0:2]=[C&H0&+0:3](-[C;H2,H1;+0:4]-[Cl,Br,I;H0;+0:6])-[O&H0&-:5]",
                "description": "tautomerization",
            },
            {
                "template": "[C;H2,H1,H0;+0:2]=[C&H0&+0:3](-[C;H2,H1;+0:4]-[Cl,Br,I;H0;+0:6])-[O&H0&-:5]>>[C;H2,H1,H0;+0:2]1-[C&H0&+0:3](-[C;H2,H1;+0:4]-1)=[O&H0&+0:5].[Cl,Br,I;H0;-:6]",
                "description": "cyclization",
            },
            {
                "template": "[C;H2,H1,H0;+0:2]1-[C&H0&+0:3](-[C;H2,H1;+0:4]-1)=[O&H0&+0:5].[R1X]>>[C;H2,H1,H0;+0:2]1-[C&H0&+0:3]([R2X])(-[C;H2,H1;+0:4]-1)-[O&H0&-1:5]",
                "description": "nucleophilic attack",
            },
            {
                "template": "[C;H2,H1,H0;+0:2]1-[C&H0&+0:3]([R2X])(-[C;H2,H1;+0:4]-1)-[O&H0&-1:5]>>[C;H2,H1,H0;-:2]-[C;H2,H1;+0:4]-[C&H0&+0:3](-[R2X])=[O&H0&+0:5]",
                "description": "hydrolysis",
            },
            {
                "template": "[C;H2,H1,H0;-:2]-[C;H2,H1;+0:4]-[C&H0&+0:3](-[R2X])=[O&H0&+0:5]>>[C;H3,H2,H1;+0:2]-[C;H2,H1;+0:4]-[C&H0&+0:3](-[R2X])=[O&H0&+0:5]",
                "description": "protonation",
            },
        ],
    },
    "Horner-Wadsworth-Emmons Olefination": {
        "scope": [
            {"name": "ester EWG", "[R1X]": "[C;H0;+0:1](=[O;H0;+0:3])(-[O;H0;+0:2])"},
            {
                "name": "acid EWG",
                "[R1X]": "[C;H0;+0:1](=[O;H0;+0:3])(-[O;H1;+0:2])",
            },
        ],
        "mechanism": [
            {
                "template": "[R1X]([C;H2;+0:4][Cl,Br,I;H0;+0:5]).[P;H0;+0:6]([O;H0;+0:7][C;H2,H1;+0:9])([O;H0;+0:12])[O;H0;+0:13]>>[R1X]([C;H2;+0:4][P;H0;+0:6](=[O;H0;+0:7])([O;H0;+0:12])[O;H0;+0:13]).[Cl,Br,I;H0;+0:5][C;H2,H1;+0:9]",
                "description": "phosphonate addition to remove halide to form active reagent",
            },
            {
                "template": "[R1X]([C;H2;+0:4][P;H0;+0:6](=[O;H0;+0:7])([O;H0;+0:12])[O;H0;+0:13])>>[R1X]([C;H1;-1:4][P;H0;+0:6](=[O;H0;+0:7])([O;H0;+0:12])[O;H0;+0:13])",
                "description": "deprotonation of phosphonate",
            },
            {
                "template": "[R1X]([C;H1;-1:4][P;H0;+0:6](=[O;H0;+0:7])([O;H0;+0:12])[O;H0;+0:13]).[C&H1&+0:10]=[O&H0&+0:11]>>[R1X](-[C&H1&+0:4](-[P&H0&+0:6](=[O&H0&+0:7])(-[O&H0&+0:12])-[O&H0&+0:13])-[C&H1&+0:10]-[O&H0&-:11])",
                "description": "conjugate addition of anion into aldehyde",
            },
            {
                "template": "[R1X](-[C&H1&+0:4](-[P&H0&+0:6](=[O&H0&+0:7])(-[O&H0&+0:12])-[O&H0&+0:13])-[C&H1&+0:10]-[O&H0&-:11])>>[R1X](-[C&H1&+0:4]1-[P&H0&+0:6](-[O&H0&-:7])(-[O&H0&+0:12])(-[O&H0&+0:13])-[O&H0&+0:11]-[C&H1&+0:10]-1)",
                "description": "formation of oxaphosphetanes",
            },
            {
                "template": "[R1X](-[C&H1&+0:4]1-[P&H0&+0:6](-[O&H0&-:7])(-[O&H0&+0:12])(-[O&H0&+0:13])-[O&H0&+0:11]-[C&H1&+0:10]-1)>>[R1X](-[C&H1&+0:4]=[C&H1&+0:10]).[P&H0&+0:6](-[O&H0&-:7])(-[O&H0&+0:12])(-[O&H0&+0:13])=[O&H0&+0:11]",
                "description": "loss of dialkyl phosphate",
            },
        ],
    },
    "Michaelis–Becker": {
        "scope": [
            {"name": "alkyl halide", "[R1X]": "[C;H2,H1;+0:2]"},
        ],
        "mechanism": [
            {
                "template": "[P;H1;+0:6](=[O;H0;+0:7])([O;H0;+0:10])[O;H0;+0:9]>>[P;H0;-1:6](=[O;H0;+0:7])([O;H0;+0:10])[O;H0;+0:9]",
                "description": "deprotonation of hydrogen phosphonate",
            },
            {
                "template": "[P;H0;-1:6](=[O;H0;+0:7])([O;H0;+0:10])[O;H0;+0:9].[Cl,Br,I;H0;+0:8][R1X]>>[P;H0;+0:6]([R1X])(=[O;H0;+0:7])([O;H0;+0:10])[O;H0;+0:9].[Cl,Br,I;H0;-1:8]",
                "description": "sn2 of phosphorous anion onto alkyl halide",
            },
        ],
    },
    "Diels-Alder": {
        "scope": [
            {"name": "ene", "[R1X]": "[C;+0:2]"},
            {"name": "imine", "[R1X]": "[N;+0:2]"},
            {"name": "oxo", "[R1X]": "[O;H0;+0:2]"},
        ],
        "mechanism": [
            {
                "template": "[C;+0:1]=[R1X].[C;+0:3]=[C;+0:4]-[C;+0:5]=[C;+0:6]>>[C&+0:1]1-[R1X]-[C;+0:3]-[C&+0:4]=[C&+0:5]-[C&+0:6]-1",
                "description": "diels-alder cycloaddition",
            },
            {
                "template": "[C;+0:1]=[R1X].[N;+0:3]=[C;+0:4]-[C;+0:5]=[C;+0:6]>>[C&+0:1]1-[R1X]-[N;+0:3]-[C&+0:4]=[C&+0:5]-[C&+0:6]-1",
                "description": "aza-diels-alder cycloaddition 1",
            },
            {
                "template": "[C;+0:1]=[R1X].[C;+0:3]=[N;+0:4]-[C;+0:5]=[C;+0:6]>>[C&+0:1]1-[R1X]-[C;+0:3]-[N&+0:4]=[C&+0:5]-[C&+0:6]-1",
                "description": "aza-diels-alder cycloaddition 2",
            },
            {
                "template": "[C;+0:3]=[N;H0;+0:4]-[C;+0:5]=[O;H0;+0:6]>>[C;+0:3]=[N;H1;+1:4]-[C;+0:5]=[O;H0;+0:6]",
                "description": "acylimine protonation",
            },
            {
                "template": "[C;+0:3]=[N;H1;+1:4]-[C;+0:5]=[O;H0;+0:6]>>[C;+1:3]-[N;H1;+0:4]-[C;+0:5]=[O;H0;+0:6]",
                "description": "acylimine tautomerization",
            },
            {
                "template": "[C;+1:3]-[N;H1;+0:4]-[C;+0:5]=[O;H0;+0:6]>>[C;+0:3]=[N;H1;+1:4]-[C;+0:5]=[O;H0;+0:6]",
                "description": "acylimine reverse tautomerization",
            },
            {
                "template": "[C;+0:1]=[R1X].[C;+0:3]=[N;H1;+1:4]-[C;+0:5]=[O;H0;+0:6]>>[C&+0:1]1-[R1X]-[C&+0:3]-[N&H1&+:4]=[C&+0:5]-[O&H0&+0:6]-1",
                "description": "acylimine cycloaddition",
            },
        ],
    },
    "Mannich-Michael": {
        "scope": [
            {
                "name": "imine p1/p2",
                "[R1X]": "[N;H1,H0;+0:2]",
                "[R2X]": "[N;H1,H0;-1:2]",
                "[R3X]": "[N;H1,H0;+0:2]",
            },
            {
                "name": "iminium p1/p2",
                "[R1X]": "[N;H2,H1;+1:2]",
                "[R2X]": "[N;H2,H1;+0:2]",
                "[R3X]": "[N;H2,H1;+1:2]",
            },
        ],
        "mechanism": [
            {
                "template": "[C;+0:1]=[R1X].[C;+0:3]=[C;+0:4]-[C;+0:5]=[C;+0:6]>>[C&+0:1](-[R2X])-[C&+0:3]-[C&+:4]-[C&+0:5]=[C&+0:6]",
                "description": "Mannich-Michael addition",
            },
            {
                "template": "[C&+0:1](-[R2X])-[C&+0:3]-[C&+:4]-[C&+0:5]=[C&+0:6]>>[C&+0:1]1-[R3X]-[C&+0:6]-[C&+0:5]=[C&+0:4]-[C&+0:3]-1",
                "description": "Mannich-Michael ring closure",
            },
        ],
    },
    "DMSO Resonance": {
        "scope": [
            {
                "name": "methyl",
                "[R1X]": "[C;H3;+0:1]",
            },
        ],
        "mechanism": [
            {
                "template": "[R1X][S;+0:2]([C:3])=[O:4]>>[R1X][S;+1:2]([O;-1:4])[C:3]",
                "description": "DMSO Resonance",
            },
            {
                "template": "[R1X][S;+1:2]([O;-1:4])[C:3].[Cl:5][P:6](=[O:7])>>[R1X][S;+0:2]([Cl:5])([O;+0:4][P:6](=[O:7]))[C:3]",
                "description": "DMSO Addition",
            },
            {
                "template": "[R1X][S;+0:2]([Cl:5])([O;+0:4][P:6](=[O:7]))[C:3]>>[R1X][S;+1:2]([Cl:5])[C:3].[O;-1:4][P:6](=[O:7])",
                "description": "DMSO chlorodimethyl sulfonium",
            },
            {
                "template": "[R1X][S;+1:2]([Cl;H0:5])[C;H3:3]>>[R1X][S;+1:2]=[C;H2:3].[Cl;H1:5]",
                "description": "chlorosulfonium tautomer",
            },
        ],
    },
    "Baeyer-Villiger": {
        "scope": [
            {
                "name": "carbon",
                "[R1X]": "[C;+0:1]",
            },
        ],
        "mechanism": [
            {
                "template": "[R1X][C:2]([C:3])=[O;H1;+1:4].[O;-1:5][P:6]=[O:7]>>[R1X][C:2]([O;+0:5][P:6]=[O:7])([C:3])[O;H1;+0:4]",
                "description": "oxidation into lactone",
            },
            {
                "template": "[R1X][C:2]([O;+0:5][P:6]=[O:7])([C:3])[O;H1;+0:4]>>[R1X][C:2]([O;+0:5][C:3])=[O;H0;+0:4].[P:6]=[O:7]",
                "description": "oxidation ring expansion",
            },
        ],
    },

    # "McMurry Coupling": {
    #     "scope": [
    #         {
    #             "name": "ketone-ketone (needs a lewis acid, 1e mechanism)",
    #             "[R1X]": "[O;H0;+0:2]",
    #             "[R2X]": "[O;H0;+0:4]",
    #         }
    #     ],
    #     "mechanism": [
    #         {
    #             "template": "[C;H0;!$(C(=O)[O;H0,H1]);+0:1]=[R1X].[C;H0;+0:3]=[R2X]>>[C;H0;!$(C(=O)[O;H0,H1]);+0:1]=[C;H0;+0:3]",
    #             "description": "McMurry ene coupling",
    #         },
    #     ],
    # },
}
