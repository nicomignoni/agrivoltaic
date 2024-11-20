using SymPy

include("model.jl")

# Angles
@syms α::real, η::real, ᾱ::real, β::real
tilt_vec = u(β) 

# Panel and crops 
@syms w::(real, positive), d::(real, positive) 
@syms c₁::real, c₂::real, c₃::real, p₁::real, p₂::real

panel = Panel(w, d, [c₁, c₂, c₃], ᾱ)
crop = Crop([p₁, p₂, 0], true) 
sun = Sun("", 0, 0, 0, α, η)

# Check the linear expression g = Au + b 
shadow = s(panel, sun, tilt_vec)
for j in 1:4
    j⁺ = 1 + (j % 4) 
    simplify(
        (g(panel, crop, sun, j, tilt_vec)) - 
        ((shadow[:,j⁺] - shadow[:,j])' * H₃' * (crop.pos - shadow[:,j]))
    ) |> println
end
