import Pkg; Pkg.activate("/home/m3g/Documentos/polyisobutileno/mddf_kb")
using Plots, LaTeXStrings, ComplexMixtures, PDBTools, EasyFit, ColorSchemes
const CM, PDB = ComplexMixtures, PDBTools

# Configurações padrão de plots
pdb_file = "p100_box.pdb"
median = "json_median"
s1 = ["3.1_100_s3_05", "3.2_100_s3_75", "3.3_100_s3_1"]
dmax_full = 3.5
plot_font = "Computer Modern"
default(fontfamily=plot_font, framestyle=:box, label=nothing, grid=true, legendfontsize=15)

# Paleta de cores seaborn_colorblind6
colors = ColorSchemes.seaborn_colorblind6.colors

# Definição de sistemas usando dicionários
systems = Dict(
    "glu_compact_05" => Dict{String,Any}(
        "pdb" => pdb_file,
        "R_glu_file" => "glu_compact_05",
        "R_water_file" => "water_compact_05"
    ),
    "glu_compact_01" => Dict{String,Any}(
        "pdb" => pdb_file,
        "R_glu_file" => "glu_compact_1",
        "R_water_file" => "water_compact_1"
    ),
    "glu_compact_75" => Dict{String,Any}(
        "pdb" => pdb_file,
        "R_glu_file" => "glu_compact_75",
        "R_water_file" => "water_compact_75"
    ),
    "glu_extended_05" => Dict{String,Any}(
        "pdb" => pdb_file,
        "R_glu_file" => "glu_extended_05",
        "R_water_file" => "water_extended_05"
    ),
    "glu_extended_01" => Dict{String,Any}(
        "pdb" => pdb_file,
        "R_glu_file" => "glu_extended_1",
        "R_water_file" => "water_extended_1"
    ),
    "glu_extended_75" => Dict{String,Any}(
        "pdb" => pdb_file,
        "R_glu_file" => "glu_extended_75",
        "R_water_file" => "water_extended_75"
    )
)

# Processar sistemas
for key in keys(systems)
    sys = systems[key]
    pdb_path = joinpath(s1[1], sys["pdb"])
    glu_file = "$median/$(sys["R_glu_file"])_average.json"
    water_file = "$median/$(sys["R_water_file"])_average.json"
    
    # Verificar existência de arquivos
    if !isfile(pdb_path) || !isfile(glu_file) || !isfile(water_file)
        println("Erro: Arquivos para $key estão faltando. Pulando...")
        continue
    end

    # Ler átomos do sistema
    sys["atoms"] = readPDB(pdb_path)

    # Carregar saída do ComplexMixtures
    sys["R_glu"] = CM.load(glu_file)
    sys["R_water"] = CM.load(water_file)

    # Calcular contribuições por resíduo da água
    rc_water = ResidueContributions(sys["R_water"], select(sys["atoms"], "resname AAMD"); dmax=dmax_full, type=:coordination_number)
    sys["rc_water"] = rc_water
    sys["rc_end_water"] = [r[end] for r in rc_water.residue_contributions]

    # Calcular contribuições por resíduo da glicose
    rc_glu = ResidueContributions(sys["R_glu"], select(sys["atoms"], "resname AAMD"); dmax=dmax_full, type=:coordination_number)
    sys["rc_glu"] = rc_glu
    sys["rc_end_glu"] = [r[end] for r in rc_glu.residue_contributions]
end

# Configuração para os gráficos
config_map = [
    ("Water/glu ext", "glu_extended", "rc_end_water", (3, 12), "A"),
    ("Water/glu comp", "glu_compact", "rc_end_water", (3, 12), "B"),
    ("Glu ext", "glu_extended", "rc_end_glu", (0, 0.8), "C"),
    ("Glu comp", "glu_compact", "rc_end_glu", (0, 0.8), "D")
]

plots = []
for (label, key, data_key, ylim, annotation) in config_map
    # Criar o gráfico base
    p = plot(systems["$(key)_01"][data_key]; xlabel="Residue", ylabel="", 
             label="$(label) 1.00 M", color=colors[1], linewidth=2, xlim=(0, 100), ylim=ylim)
    plot!(systems["$(key)_75"][data_key]; label="$(label) 0.75 M", color=colors[2], linewidth=2)
    plot!(systems["$(key)_05"][data_key]; label="$(label) 0.50 M", color=colors[3], linewidth=2, legend=:topright)

    # Calcular deslocamento ajustado para posicionar a anotação mais alta
    deslocamento = (ylim[2] - ylim[1]) * 0.9  # 90% do intervalo para subir a anotação
    posicao_y = ylim[1] + deslocamento       # Nova posição Y da anotação

    # Adicionar anotação no gráfico
    annotate!([(5.0, posicao_y, text(annotation, plot_font, 30, :left, color="#000000"))])

    # Adicionar o gráfico final à lista
    push!(plots, p)
end

# Combinar gráficos em um layout de 2x2
layout = @layout [a b; c d]
final_plot = plot(plots..., layout=layout, size=(1200, 800),
                  ylabel=L"\mathrm{Coordination~number~up~to~3.5~Å}",
                  xtickfontsize=25, ytickfontsize=25, xlabelfontsize=20, ylabelfontsize=15,
                  margin=8Plots.Measures.mm, dpi=500, legendfontsize=15, legend=:topright)

# Salvar o gráfico final
savefig(final_plot, "imagens/coordenation_layout_2x2_$dmax_full.png")
