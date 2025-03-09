source('bases on [l,u] with tuned gb.R')

integrate(function(t) S1(t)^2*gb(t), l, u)$value |> round(10)
integrate(function(t) T1_normed(t)^2*gb(t), l, u)$value |> round(10)
integrate(function(t) T2_normed(t)^2*gb(t), l, u)$value |> round(10)
integrate(function(t) T3_normed(t)^2*gb(t), l, u)$value |> round(10)
integrate(function(t) T4_normed(t)^2*gb(t), l, u)$value |> round(10)
integrate(function(t) T5_normed(t)^2*gb(t), l, u)$value |> round(10)

integrate(function(t) S1(t)*gb(t), l, u)$value |> round(10)
integrate(function(t) T1_normed(t)*gb(t), l, u)$value |> round(10)
integrate(function(t) T2_normed(t)*gb(t), l, u)$value |> round(10)
integrate(function(t) T3_normed(t)*gb(t), l, u)$value |> round(10)
integrate(function(t) T4_normed(t)*gb(t), l, u)$value |> round(10)
integrate(function(t) T5_normed(t)*gb(t), l, u)$value |> round(10)

integrate(function(t) T1_normed(t)*S1(t)*gb(t), l, u)$value |> round(10)
integrate(function(t) T2_normed(t)*S1(t)*gb(t), l, u)$value |> round(10)
integrate(function(t) T3_normed(t)*S1(t)*gb(t), l, u)$value |> round(10)
integrate(function(t) T4_normed(t)*S1(t)*gb(t), l, u)$value |> round(10)
integrate(function(t) T5_normed(t)*S1(t)*gb(t), l, u)$value |> round(10)

integrate(function(t) T2_normed(t)*T1_normed(t)*gb(t), l, u)$value |> round(10)
integrate(function(t) T3_normed(t)*T1_normed(t)*gb(t), l, u)$value |> round(10)
integrate(function(t) T4_normed(t)*T1_normed(t)*gb(t), l, u)$value |> round(10)
integrate(function(t) T5_normed(t)*T1_normed(t)*gb(t), l, u)$value |> round(10)

integrate(function(t) T3_normed(t)*T2_normed(t)*gb(t), l, u)$value |> round(10)
integrate(function(t) T4_normed(t)*T2_normed(t)*gb(t), l, u)$value |> round(10)
integrate(function(t) T5_normed(t)*T2_normed(t)*gb(t), l, u)$value |> round(10)

integrate(function(t) T4_normed(t)*T3_normed(t)*gb(t), l, u)$value |> round(10)
integrate(function(t) T5_normed(t)*T3_normed(t)*gb(t), l, u)$value |> round(10)

integrate(function(t) T5_normed(t)*T4_normed(t)*gb(t), l, u)$value |> round(10)

