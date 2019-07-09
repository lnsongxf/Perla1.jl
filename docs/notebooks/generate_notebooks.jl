using Pkg
using Weave

pkg"activate ../.." # activate main Perla1 repository

fileset = readdir(@__DIR__)
for file in fileset
    if occursin(".jmd", file) 
        Weave.notebook(file, :pwd, -1, "--allow-errors") # executes also
        println("Weaved $file")
    end
end
