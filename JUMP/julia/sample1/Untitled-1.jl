using JuMP
import Ipopt
using Plots

rocket = Model(Ipopt.Optimizer)
set_silent(rocket)

h_0 = 1    # Initial height
v_0 = 0    # Initial velocity
m_0 = 1    # Initial mass
g_0 = 1    # Gravity at the surface

T_c = 3.5  # Used for thrust
h_c = 500  # Used for drag
v_c = 620  # Used for drag
m_c = 0.6  # Fraction of initial mass left at end

c = 0.5 * sqrt(g_0 * h_0)  # Thrust-to-fuel mass
m_f = m_c * m_0              # Final mass
D_c = 0.5 * v_c * m_0 / g_0  # Drag scaling
T_max = T_c * g_0 * m_0        # Maximum thrust

n = 800    # Time steps

@variables(rocket, begin
    Δt ≥ 0, (start = 1 / n) # Time step
    # State variables
    v[1:n] ≥ 0            # Velocity
    h[1:n] ≥ h_0          # Height
    m_f ≤ m[1:n] ≤ m_0    # Mass
    # Control variables
    0 ≤ T[1:n] ≤ T_max    # Thrust
end);

@objective(rocket, Max, h[n])

fix(v[1], v_0; force = true)
fix(h[1], h_0; force = true)
fix(m[1], m_0; force = true)
fix(m[n], m_f; force = true)

@NLexpressions(
    rocket,
    begin
        # Drag(h,v) = Dc v^2 exp( -hc * (h - h0) / h0 )
        drag[j = 1:n], D_c * (v[j]^2) * exp(-h_c * (h[j] - h_0) / h_0)
        # Grav(h)   = go * (h0 / h)^2
        grav[j = 1:n], g_0 * (h_0 / h[j])^2
        # Time of flight
        t_f, Δt * n
    end
);

for j in 2:n
    # h' = v
    # Rectangular integration
    # @NLconstraint(rocket, h[j] == h[j - 1] + Δt * v[j - 1])
    # Trapezoidal integration
    @NLconstraint(rocket, h[j] == h[j-1] + 0.5 * Δt * (v[j] + v[j-1]))
    # v' = (T-D(h,v))/m - g(h)
    # Rectangular integration
    # @NLconstraint(
    #     rocket,
    #     v[j] == v[j - 1] + Δt *((T[j - 1] - drag[j - 1]) / m[j - 1] - grav[j - 1])
    # )
    # Trapezoidal integration
    @NLconstraint(
        rocket,
        v[j] ==
        v[j-1] +
        0.5 *
        Δt *
        (
            (T[j] - drag[j] - m[j] * grav[j]) / m[j] +
            (T[j-1] - drag[j-1] - m[j-1] * grav[j-1]) / m[j-1]
        )
    )
    # m' = -T/c
    # Rectangular integration
    # @NLconstraint(rocket, m[j] == m[j - 1] - Δt * T[j - 1] / c)
    # Trapezoidal integration
    @NLconstraint(rocket, m[j] == m[j-1] - 0.5 * Δt * (T[j] + T[j-1]) / c)
end

println("Solving...")
optimize!(rocket)
solution_summary(rocket)

println("Max height: ", objective_value(rocket))
function my_plot(y, ylabel)
    return Plots.plot(
        (1:n) * value.(Δt),
        value.(y)[:];
        xlabel = "Time (s)",
        ylabel = ylabel,
    )
end

Plots.plot(
    my_plot(h, "Altitude"),
    my_plot(m, "Mass"),
    my_plot(v, "Velocity"),
    my_plot(T, "Thrust");
    layout = (2, 2),
    legend = false,
    margin = 1Plots.cm,
)