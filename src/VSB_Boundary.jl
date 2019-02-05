"""
`Boundary`


"""
mutable struct Boundary

    bodyPTS::NUMB_MATRIX
    tHats::NUMB_MATRIX
    nHats::NUMB_MATRIX
    NPTS_BODY::Int

    function Boundary(body_pts::NUMB_MATRIX,
                      t_hats::NUMB_MATRIX,
                      n_hats::NUMB_MATRIX)
        new(body_pts,t_hats,n_hats,length(body))
    end
end
