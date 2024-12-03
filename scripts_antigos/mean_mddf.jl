import Pkg; Pkg.activate("/home/m3g/Documentos/1_analysis/1_analises/p100_analysis")
using Plots, LaTeXStrings, ComplexMixtures, PDBTools, EasyFit, ColorSchemes
const CM, PDB = ComplexMixtures, PDBTools

# Diretórios principais
main_dirs = ["glu_compact", "glu_extended", "water_compact", "water_extended"]
sub_dirs = ["1", "05", "75"]

# Função para carregar e plotar mddf por resíduo
function plot_mddf_by_residue(dir, subdir)
    # Construindo o caminho para o diretório específico
    dir_path = joinpath(dir, subdir)
    
    # Encontrar o arquivo JSON e o arquivo PDB no subdiretório
    json_files = filter(f -> endswith(f, ".json"), readdir(dir_path, join=true))
    pdb_files = filter(f -> endswith(f, ".pdb"), readdir(dir_path, join=true))
    
    # Verifica se os arquivos existem e se há apenas um de cada tipo
    if length(json_files) == 1 && length(pdb_files) == 1
        file_json = json_files[1]
        file_pdb = pdb_files[1]
        
        # Carrega os arquivos JSON e PDB
        data = CM.load(file_json)
        atoms = PDBTools.readPDB(file_pdb)
        
        # Seleciona as partes relevantes do sistema
        protein = PDBTools.select(atoms, "resname AAMD")
        solute = CM.AtomSelection(protein, nmols=1)
        
        bgl = PDBTools.select(atoms, by = at -> PDBTools.resname(at) in ["AGL","BGL"])
        solvent = CM.AtomSelection(bgl, natomspermol=24)
        
        # Seleção de resíduos específicos
        C = PDBTools.select(bgl, by = at -> PDBTools.element(at) == "C")
        data_C_mddf = CM.contributions(solvent, data.solvent_atom, C)
        data_C_mddf = movingaverage(data_C_mddf, 10).x
        
        H = PDBTools.select(bgl, by = at -> PDBTools.element(at) == "H")
        data_H_mddf = CM.contributions(solvent, data.solvent_atom, H)
        data_H_mddf = movingaverage(data_H_mddf, 10).x
        
        O = PDBTools.select(bgl, by = at -> PDBTools.element(at) == "O")
        data_O_mddf = CM.contributions(solvent, data.solvent_atom, O)
        data_O_mddf = movingaverage(data_O_mddf, 10).x
        
        # Plotando os dados
        x = data.d
        y_total = data.mddf
        
        plot(x, y_total, framestyle = :box, label = "Total", line = :solid, linewidth = 4, linecolor = :blue)
        plot!(x, data_C_mddf, label = "C", line = :solid, linewidth = 4, linecolor = :red)
        plot!(x, data_H_mddf, label = "H", line = :solid, linewidth = 4, linecolor = :green)
        plot!(x, data_O_mddf, label = "O", line = :solid, linewidth = 4, linecolor = :orange)
        
        # Configurações finais do plot
        plot!(title = "MDDF by Residue", xlabel = "Distance (Å)", ylabel = "g(r)", xlim = [1.0, 10.0], ylim = [0.0, 5.0])
        
        # Salvando o gráfico
        output_file = joinpath(dir_path, "mddf_by_residue_plot.png")
        savefig(output_file)
    else
        println("Arquivos não encontrados ou múltiplos arquivos em $dir/$subdir")
    end
end

# Itera sobre os diretórios principais e subdiretórios
for dir in main_dirs
    for subdir in sub_dirs
        plot_mddf_by_residue(dir, subdir)
    end
end
