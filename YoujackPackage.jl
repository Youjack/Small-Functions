module YoujackPackage

using Plots, DataFrames, GLM

export linearfitplot

youjack_plotcolor = let color=palette(:tab10), len=length(color)
    reshape([color[i] for i in 1:len], (1,len))
end

"""
    linearfitplot(
        labels::Array, data::DataFrame,
        xrange::Tuple, axeslabel::Tuple;
        basis = true)
"""
function linearfitplot(
    labels::Array{String,1}, data::DataFrame,
    xrange::NTuple{2,Real}, axeslabel::NTuple{2,String};
    basis::Bool = true
)
    len = length(labels)

    # fit
    names = propertynames(data)
    coefs = DataFrame(intercept=Float64[], slope=Float64[], R²=Float64[])
    lines = Array{Function,1}(undef,len)
    for i in 1:len
        model = lm(
            FormulaTerm(Term(names[2i]), (ConstantTerm(basis ? 1 : 0), Term(names[2i-1])) ),
            data[!,2i-1:2i] )
        push!(coefs, (
            basis ? coef(model)[1] : 0, 
            basis ? coef(model)[2] : coef(model)[1], 
            r2(model)) )
        lines[i] = ( x -> coefs[i,1] + coefs[i,2] * x )
    end
    
    # plot
    gr()
    fig = plot(legend=:outerright, framestyle=:origin,
        xlims=xrange,
        xlabel=axeslabel[1],ylabel=axeslabel[2])
    for i in 1:len
        scatter!(fig, data[!,2i-1],data[!,2i], 
            label=labels[i], seriescolor=youjack_plotcolor[i])
        plot!(fig, lines[i], xrange[1],xrange[2], 
            label=false, seriescolor=youjack_plotcolor[i])
    end
    
    # return
    (coefs,fig)
end

end # module YoujackPackage