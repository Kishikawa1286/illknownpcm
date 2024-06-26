# 修正不要なら a を返す
# a と b が十分に近い場合に b を返す
@inline function correctPrecisionLoss(a, b)
    if abs(a - b) < 1e-8
        return b
    end
    return a
end

@inline function nearlyEqual(a::Number, b::Number)
    if abs(a - b) < 1e-8
        return true
    end
    return false
end

@inline function nearlyEqualLoose(a::Number, b::Number)
    if abs(a - b) < 1e-3
        return true
    end
    return false
end

