module TextIO_Utils
    using Printf
    import LinearAlgebra: eigen


    function Write_tan_res(Output_file, x, W)
        write(Output_file, "x: \n")
        Write_pretty_vector(Output_file, x)
        write(Output_file, "W: \n")
        Write_pretty_matrix(Output_file, W)
        write(Output_file, "############### \n")
        eigs = eigen(W)
        write(Output_file, "values: \n ")
        Write_pretty_vector(Output_file, eigs.values)
        write(Output_file, "vectors: \n")
        Write_pretty_matrix(Output_file, eigs.vectors)
    end

    function Write_pretty_matrix(io, mat)
        for i in 1:size(mat)[1]
            for j in 1:size(mat)[2]
                if typeof(mat[i, j]) == ComplexF64
                    write(io, @sprintf("%20.15s", @sprintf("%5.3g", real(mat[i, j])) * " + " * @sprintf("%5.3g", imag(mat[i, j])) * "i  "))
                else
                    write(io, @sprintf("%10.3g  ", mat[i, j]))
                end
            end
            write(io, "\n")
        end
    end

    function Write_pretty_vector(io, vec)
        for i in 1:size(vec)[1]
            if typeof(vec[i]) == ComplexF64
                write(io, @sprintf("%20.15s \n", @sprintf("%5.3g", real(vec[i])) * " + " * @sprintf("%5.3g", imag(vec[i])) * "i "))
            else
                write(io, @sprintf("%10.5g  \n", vec[i]))
            end
            write(io, "\n")
        end
    end
    
    function Print_pretty_matrix(mat)
        for i in 1:size(mat)[1]
            for j in 1:size(mat)[2]
                print(@sprintf("%10.5g  ", mat[i, j]))
            end
            print("\n")
        end
    end
end