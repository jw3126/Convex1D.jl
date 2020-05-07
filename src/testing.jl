module Testing
export testdatapath

function testdatapath(args...)
    dir = normpath(joinpath(@__DIR__, "..", "test", "data"))
    joinpath(dir, args...)
end

end#module
