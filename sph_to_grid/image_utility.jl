@inline function check_in_image(pos::Float64, hsml::Float64,
                                minCoord::Float64, maxCoord::Float64)

    if ( (minCoord - hsml) <= pos <= (maxCoord + hsml) )
        return true
    else
        return false
    end

end

@inline function find_min_pixel(pos::Float64, hsml::Float64,
                                minCoord::Float64, pixSize::Float64)

    pixmin = floor((pos - hsml - minCoord) / pixSize )
    if pixmin < 1
        return 1
    else
        return pixmin
    end
end

@inline function find_max_pixel(pos::Float64, hsml::Float64,
                                minCoord::Float64, pixSize::Float64,
                                max::Int64)

    pixmax = floor((pos + hsml - minCoord) / pixSize )
    if pixmax > max
        return max
    else
        return pixmax
    end
end
