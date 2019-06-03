using Pkg 
pkg"add https://github.com/arnavs/Weave.jl" # version of Weave that supports notebook
using Weave

pkg"activate .." # activate main Perla1 repository

fileset = readdir(@__DIR__)
for file in fileset
    if occursin(".jmd", file)
        Weave.notebook(file) # executes also
        println("Weaved $file")
    end
end
