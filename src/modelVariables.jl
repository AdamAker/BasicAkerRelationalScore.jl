using BasicAkerRelationalScore
using OrderedCollections
using LaTeXStrings

function splitDataFrame(bigDataFrame,featureNames)
    
    featuresDataFrame = bigDataFrame[:,featureNames]
    targetsDataFrame = bigDataFrame[:,Not(featureNames)]

    return featuresDataFrame,targetsDataFrame

end

function makeFeaturesDict(featuresDataFrame,targetsDataFrame,acceptance)

    featuresDict = OrderedDict()
    featureList = []
    powerBARSs = []

    for featureName∈propertynames(featuresDataFrame)

        df = DataFrame()
        df[:,featureName]=featuresDataFrame[:,featureName]
        df = hcat(df,targetsDataFrame)
        
        targetsDict = BasicAkerRelationalScore.makeTargetsDict(df,featureName,acceptance)
        targetsDict = BasicAkerRelationalScore.calcPBARS(targetsDict)

        featuresDict[featureName] = targetsDict
        push!(powerBARSs,targetsDict[:pBARS])
        push!(featureList, featureName)

    end

    featuresDict[:powerBARSs] = powerBARSs
    featuresDict[:featureList] = propertynames(featuresDataFrame)

    return featuresDict

end

function makeFeaturesDict(dataFrame,acceptance)

    α=acceptance
    n=length(propertynames(dataFrame))
    BARSMatrix = zeros(n,n)
    R²Matrix = zeros(n,n)
    featuresDict = OrderedDict()
    featureList = []
    powerBARSs = []
    for i∈1:n
        featureName = propertynames(dataFrame)[i]
        push!(featureList, featureName)
        for j∈1:n
            targetName = propertynames(dataFrame)[j]
            if !isequal(featureName,targetName)
                targetsDict = BasicAkerRelationalScore.makeTargetsDict(dataFrame,featureName,α)
                featuresDict[featureName] = targetsDict
                BARSMatrix[n+1-i,j] = targetsDict[targetName][:BARS]
                R²Matrix[n+1-i,j] = targetsDict[targetName][:R²]
            else
                BARSMatrix[n+1-i,j]= 1.0
                R²Matrix[n+1-i,j]=1.0
            end
        end

    end

    pMatrix = R²Matrix.*BARSMatrix

    featuresDict[:BARSMatrix] = BARSMatrix
    featuresDict[:R²Matrix] = R²Matrix
    featuresDict[:pMatrix] = pMatrix

    return featuresDict

end


function plotPowerBARSs(featuresDict)

    finalIndex = length(featuresDict[:powerBARSs])
    sortedPowerBARSs = [0.0 for i∈1:finalIndex]
    sortedFeatureNames = ["0" for i∈1:finalIndex]
    sortedIndicies = sortperm(featuresDict[:powerBARSs])

    for sortingIndex ∈ sortedIndicies
        sortedPowerBARSs[finalIndex] = featuresDict[:powerBARSs][sortingIndex]
        sortedFeatureNames[finalIndex] = string(featuresDict[:featureList][sortingIndex])
        finalIndex-=1
    end

    PowerBARSsPlot = Plots.bar(sortedFeatureNames,
                        sortedPowerBARSs,
                        title = "PBARSs of the Features",
                        label = L"PBARS=\sum_{i=1}^{N}R_{i}^2BARS_i",
                        xticks = :all,
                        xrotation = 60.0)
    xlabel!("Feature Name")
    ylabel!("PowerBARS")

    return PowerBARSsPlot

end

function generateModelvars(featuresDict)

    D = OrderedDict(zip(featuresDict[:featureList],featuresDict[:powerBARSs]))
    sortedD = sort(D, byvalue = true, rev=true)

    modelVariables = OrderedDict()
    for keyNames∈ sortedD.keys
        if length(featuresDict[keyNames][:targetsAbove])>0
            modelVariables[keyNames] = featuresDict[keyNames][:targetsAbove]
        end
    end

    return modelVariables

end

function matrixPlot(name,matrix,dataFrame,cutoff,color,cScheme,fontSize)
    n = length(propertynames(dataFrame))

    hmap=plot(heatmap(
        xticks = (1:1:n, propertynames(dataFrame)),
        yticks = (n:-1:1, propertynames(dataFrame)),
        matrix, 
        title = name*" Heat Map",
        xlabel = "targets",
        ylabel = "features",
        xrotation = 60.0,
        c=cScheme
    ))

    for tick∈1:n
        vline!([tick+.5], color=color, label=false)
        hline!([tick+.5], color=color, label=false)
    end

    xpoints = [ i for i∈1:n]
    ypoints = [ j for j∈1:n]
    for ycoor∈ypoints
        for xcoor∈xpoints
            value = matrix[ycoor,xcoor]
            if value<cutoff
                value = 0
            end
            plot!((xcoor,ycoor),
                legend = false,
                series_annotation=text.(round((value),sigdigits=2),fontSize,color)
            )
        end
    end

    return hmap
end

