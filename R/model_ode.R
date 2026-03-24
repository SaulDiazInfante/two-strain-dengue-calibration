#' Internal ODE System for the Two-Strain Dengue Model
#'
#' @param time Simulation time supplied by [deSolve::ode()].
#' @param state Named numeric state vector.
#' @param parameters Named numeric parameter vector.
#'
#' @return A list containing the derivatives expected by [deSolve::ode()].
#' @keywords internal
#' @noRd
dengue_model_ode <- function(time, state, parameters) {
    with(as.list(c(state, parameters)), {
        total_humans <- S + I1 + I2 + S1 + Y1h + Y1c + R

        force_mosquito_i1 <- (betaM * b / total_humans) * I1
        force_mosquito_i2 <- (betaM * b / total_humans) * I2
        force_mosquito_y1h <- (betaM * b / total_humans) * Y1h
        force_mosquito_y1c <- (betaM * b / total_humans) * Y1c
        total_mosquito_force <- force_mosquito_i1 + force_mosquito_i2 +
            force_mosquito_y1h + force_mosquito_y1c

        force_human_m1 <- (betaH * b / total_humans) * M1
        force_human_m2 <- (betaH * b / total_humans) * M2

        dMs <- LambdaM - total_mosquito_force * Ms - muM * Ms
        dM1 <- force_mosquito_i1 * Ms - muM * M1
        dM2 <- (
            force_mosquito_i2 +
                force_mosquito_y1h +
                force_mosquito_y1c
        ) * Ms - muM * M2
        dS <- LambdaS - (force_human_m1 + force_human_m2) * S - muH * S
        dI1 <- force_human_m1 * S - (alphaC + muH) * I1
        dI2 <- force_human_m2 * S - (alphaC + muH) * I2
        dS1 <- Lambda1 - vsigma * force_human_m2 * S1 - muH * S1
        dY1c <- (1 - vtheta) * vsigma * force_human_m2 * S1 - (alphaC + muH) * Y1c
        dY1h <- vtheta * vsigma * force_human_m2 * S1 - (alphaH + muH) * Y1h
        dR <- alphaC * (I1 + I2 + Y1c) + alphaH * Y1h - muH * R
        dz <- p * (dI1 + dI2 + dY1c)

        list(c(dMs, dM1, dM2, dS, dI1, dI2, dS1, dY1c, dY1h, dR, dz))
    })
}

FuncionODE <- dengue_model_ode
