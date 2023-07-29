``` julia
import Pkg
using Plots, LaTeXStrings, ComplexMixtures, PDBTools, EasyFit, ColorSchemes
const CM=ComplexMixtures
const PDB=PDBTools

# simulacao acrem glicose na concentracao x.xx

# replicas 
files_water_pure  =  [ load("./ws1_w.json"), load("./ws2_w.json"), load("./ws3_w.json") ]    # agua pura
files_water_glu   =  [ load("./ws1_05.json"), load("./ws2_05.json"), load("./ws3_05.json") ] # agua + glicose
files_glucose     =  [ load("./gs1_05.json"), load("./gs2_05.json"), load("./gs3_05.json") ] # glicose

# criando o vetor
vec_water_pure    = []
vec_water_glu     = []
vec_glucose       = []

# juntando arquivos 
vec_water_pure    = merge(files_water_rep)
vec_water_glu     = merge(files_water_p100)
vec_glucose       = merge(files_glucose)

# salvando os resultado em um único json
save(vec_water_pure,"water_full.json")
save(vec_water_glu,"water_glu_05.json")
save(vec_glucose,"glucose_full_05.json")

function read_data()

    # carregando os arquivos "json"
    data1 = CM.load("water_full.json")      # agua
    data2 = CM.load("water_glu_05.json")    # agua + glucose
    data3 = CM.load("glucose_full_05.json") # glucose + agua
    data4 = CM.load("gs1_05.json")          # s1 - 0.50
    data5 = CM.load("gs2_05.json")          # s2 - 0.50
    data6 = CM.load("gs3_05.json")          # s3 - 0.50
    
    return data1, data2, data3, data4, data5, data6
end

function myplots(data1,data2,data3,data4,data5,data6)
    
    # configuracoes padrao de plots
    
    plot_font = "Computer Modern"
    default(fontfamily=plot_font, framestyle=:box, label=nothing, grid=true, legendfontsize=12)

    # aplicando a média móvel nos dados de mddf para os arquivos json completos (s1+s2+s3)
    
    data1_mddf = movingaverage(data1.mddf,10).x  # mddf - agua pura
    data2_mddf = movingaverage(data2.mddf,10).x  # mddf - agua pura + glicose
    data3_mddf = movingaverage(data3.mddf,10).x  # mddg - glicose
    
    # -------------------------------------------------------
    
    # plotando mddf
    
    x = data1.d
    y = data1_mddf

    plot((x),(y),
            framestyle      =:box,
            title           ="",
            label           =L"\textrm{pure~water}",
            line            =:solid,
            ticks           =true,
            linewidth       =3,
            linecolor       =RGB(1.0,0.41,0.71), #hotpink
            xaxis           =true,
            yaxis           =true
    )
    
    x = data2.d
    y = data2_mddf
    
    plot!((x),(y),
            framestyle      =:box,
            title           ="",
            label           =L"\textrm{water}",
            line            =:solid,
            linealpha       =0.7,
            ticks           =true,
            linewidth       =3,
            linecolor       =RGB(0.29, 0.0, 0.51), # indigo
            xaxis           =true, 
            yaxis           =true
    )
    
    x = data3.d
    y = data3_mddf
    
    plot!((x),(y),
            framestyle      =:box,
            title           = "",
            label           = L"\mathrm{glucose~1.0~mol~L}^{-1}",
            line            = :solid,
            ticks           = true,
            linewidth       = 3,
            linecolor       = RGB(0.93, 0.46, 0.13), # blue2
            xaxis           = true, 
            yaxis           = true
    )
    
    mddf_all = plot!(title  = "",
            xaxis           = ("r/Å",0:1.5:15),
            ylabel          = L"\textrm{Minimum~distance~distribution~functions}",
            dpi             = 1000
    )
	
	
	# salvando dados da mddf
	
	# -------------------------------------------------------
	savefig(mddf_all,"./mddf_100_05.png")
	# -------------------------------------------------------
	
	# mddf por atomo
	
     atoms = PDB.readPDB("p100_box.pdb")
    
     protein = PDB.select(atoms,"resname AAMD")
     solute = CM.Selection(protein,nmols=1)
    
     bgl = PDB.select(atoms, by = at -> resname(at) in ["AGL","BGL"])
     solvent = CM.Selection(bgl,natomspermol=24)
    
     # Configurando as seleções 
    
     C = PDB.select(bgl, by = at -> element(at) == "C")
     data3_C_mddf = contributions(solvent,data3.solvent_atom,C) 
     data3_C_mddf = movingaverage(data3_C_mddf,10).x
    
     H = PDB.select(bgl, by = at -> element(at) == "H")
     data3_H_mddf = contributions(solvent,data3.solvent_atom,H)
     data3_H_mddf = movingaverage(data3_H_mddf,10).x
        
     O = PDB.select(bgl, by = at -> element(at) == "O")
     data3_O_mddf = contributions(solvent,data3.solvent_atom,O)
     data3_O_mddf = movingaverage(data3_O_mddf,10).x
    
     # -------------------------------------------------------
    
     x = data3.d
     y = data3_mddf
    
     plot((x),(y),
            framestyle     =:box,
            title          = "",
            label          = L"\textrm{total}",
            line           = :solid,
            ticks          = true,
            linewidth      = 3,
            linecolor      = RGB(1.0,0.41,0.71), #hotpink
            xaxis          = true, 
            yaxis          = true
    )
    
    x = data3.d
    y = data3_H_mddf
    
    plot!((x),(y),
            framestyle     =:box,
            title          = "",
            label          = L"\textrm{H}",
            line           = :solid,
            ticks          = true,
            linewidth      = 3,
            linecolor      = RGB(0.29, 0.0, 0.51), # indigo
            xaxis          = true, 
            yaxis          = true
    )
    
    x = data3.d
    y = data3_O_mddf
    
    plot!((x),(y),
            framestyle     =:box,
            title          = "",
            label          = L"\textrm{O}",
            line           = :solid,
            grid           = true,
            ticks          = true,
            linewidth      = 3,
            linecolor      = RGB(0.93, 0.46, 0.13), # blue2
            xaxis          = true, 
            yaxis          = true
    )
    
    mddf_atoms_glu = plot!(title    = "",
            titlefontsize           = 12,
            xaxis                   = ("r/Å",0:1.5:15),
            ylabel                  = L"\textrm{Minimum-distance~distribution~functions}",
            dpi                     = 1000
    )

	# salvando dados da mddf por atomo

	# -------------------------------------------------------
	savefig(mddf_atoms_glu,"./mddf_atoms.png")
	# -------------------------------------------------------

	# ESTUDO DAS INTEGRAIS DE KIRKWOOD-BUFF
	
     data1_kb = data1.kb  # agua - full
     data2_kb = data2.kb  # agua + glicose - full
     data3_kb = data3.kb  # glicose - full
     data4_kb = data4.kb  # s1 - glicose
     data5_kb = data5.kb  # s2 - glicose
     data6_kb = data6.kb  # s3 - glicose
    
     data1_kb = movingaverage(data1_kb,10).x
     data2_kb = movingaverage(data2_kb,10).x
     data3_kb = movingaverage(data3_kb,10).x
     data4_kb = movingaverage(data4_kb,10).x
     data5_kb = movingaverage(data5_kb,10).x
     data6_kb = movingaverage(data6_kb,10).x
    
     x = data2.d
     y = data2_kb 
    
     plot((x),(y)/1000,
            framestyle       =:box,
            title            = "",
            label            = L"\mathrm{water}",
            line             = :solid,
            linealpha        = 0.9,
            ticks            = true,
            linewidth        = 3,
            linecolor        = RGB(0.29, 0.0, 0.51), # indigo
            xaxis            = true, 
            yaxis            = true
     )
    
     x = data3.d
     y = data3_kb 
    
     plot!((x),(y)/1000,
            framestyle      =:box,
            title           = "",
            label           = L"\mathrm{glucose~0.5~mol~L}^{-1}",
            line            = :solid,
            linealpha       = 0.9,
            ticks           = true,
            linewidth       = 3,
            linecolor       = RGB(1.0,0.41,0.71), #hotpink
            xaxis           = true, 
            yaxis           = true
     )
    
     kb_water_glu            = plot!(title = "",
    	       legend          =:bottomright,
            xaxis           = (L"\textrm{r/Å}",0:1.5:15),
            ylabel          = L"\mathrm{Kirkwood~Buff}",
            dpi             = 1000
     )
    
     # salvando dados das integrais de KB
    
     # -------------------------------------------------------
     savefig(kb_water_glu,"./kb_1.png")
     # -------------------------------------------------------

     # plotando todas as concentrações

	x = data2.d
	y = data2_kb

	plot((x),(y)/1000,
            framestyle       =:box,
            title            = "",
            label            = L"\mathrm{water}",
            line             = :solid,
            linealpha        = 0.9,
            ticks            = true,
            linewidth        = 3,
            linecolor        = RGB(1.0,0.41,0.71), #hotpink
            xaxis            = true, 
            yaxis            = true
    )

	x = data3.d
	y = data3_kb

	plot!((x),(y)/1000,
            framestyle      =:box,
            title           = "",
            label           = L"\mathrm{glucose~0.5~mol~L}^{-1}",
            line            = :solid,
            linealpha       = 0.9,
            ticks           = true,
            linewidth       = 3,
            linecolor       = :red,
            xaxis           = true, 
            yaxis           = true
    )

	x = data4.d
	y = data4_kb

	plot!((x),(y)/1000,
            framestyle       =:box,
            title            = "",
            label            = L"\mathrm{Polymer~minimized}",
            line             = :solid,
            linealpha        = 0.9,
            ticks            = true,
            linewidth        = 3,
            linecolor        = RGB(1.0, 0.73, 0.06), # dark goldenrod1
            xaxis            = true, 
            yaxis            = true
    )

	x = data5.d
	y = data5_kb

	plot!((x),(y)/1000,
            framestyle       =:box,
            title            = "",
            label            = L"\mathrm{Polymer~open}",
            line             = :solid,
            linealpha        = 0.9,
            ticks            = true,
            linewidth        = 3,
            linecolor        = RGB(0.0, 1.0, 0.0), # green1
            xaxis            = true, 
            yaxis            = true
    )

	x = data6.d
	y = data6_kb

	plot!((x),(y)/1000,
            framestyle      =:box,
            title           = "",
            label           = L"\mathrm{Polymer~collapsed}",
            line            = :solid,
            linealpha       = 0.9,
            ticks           = true,
            linewidth       = 3,
            linecolor       = RGB(0.55,0.0,0.55), # magenta4
            xaxis           = true, 
            yaxis           = true
    )

	kb_mddf_all_05 = plot!(title  = L"\textrm{All replicas}",
            xaxis                  = ("r/Å",0:1.5:15),
            ylabel                 = L"\textrm{Kirkwood~Buff}",
            legend                 =(:bottomright),
            dpi                    = 500
	)

	# salvando dados das integrais de KB

	# -------------------------------------------------------
	savefig(kb_mddf_all_05,"./kb_100_all_05.png")
	# -------------------------------------------------------

	return mddf_all, mddf_atoms_glu, kb_water_glu, kb_mddf_all_05

end

```
