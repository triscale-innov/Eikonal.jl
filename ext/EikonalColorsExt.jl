module EikonalColorsExt

using Eikonal
using Colors

function Eikonal.img2array(img, color_mapping)
    dict = map(color_mapping) do (color, value)
        parse(Colorant, color) => value
    end |> Dict

    col = collect(keys(dict))
    map(img) do c
        (_, i) = findmin(c′->colordiff(c,c′), col)
        dict[col[i]]
    end
end

end
