from flask import Flask, jsonify, send_from_directory, request
import os
try:
    import google.generativeai as genai
    _gemini_available = True
except ImportError:
    _gemini_available = False


app = Flask(__name__)

# Data storage
data = {
    "11": {
        "Physics": {
            "Physics & Measurement": {
                "difficulty": "Easy",
                "description": "Physics & Measurement is the gateway chapter of JEE Physics. It establishes the language of physics — units, dimensions, and the art of precise measurement. Every physical quantity in nature can be expressed in terms of seven fundamental SI units. Dimensional analysis is a powerful tool to derive formulas, check equations, and convert units. Understanding significant figures and errors is essential for experimental accuracy, which is tested extensively in JEE practical-based questions.",
                "roadmap": ["Physical Quantities & SI Units", "Dimensional Analysis", "Significant Figures", "Errors in Measurement", "Instruments — Vernier Caliper & Screw Gauge"],
                "topics": ["Fundamental & Derived Quantities", "Dimensions of Physical Quantities", "Principle of Homogeneity", "Absolute, Relative & Percentage Error", "Propagation of Errors", "Least Count of Instruments", "Order of Magnitude"],
                "formulas": ["[M L T] notation for dimensions", "% Error = (ΔA/A) × 100", "If Z = A×B, then ΔZ/Z = ΔA/A + ΔB/B", "Least count = MSD − VSD"],
                "laws": ["Principle of Homogeneity — dimensions on both sides of an equation must be equal"],
                "applications": ["Checking correctness of physical equations", "Deriving unknown formulas", "Converting units between systems"],
                "common_mistakes": ["Forgetting to square dimensions in area/volume", "Adding quantities with different dimensions", "Confusing precision with accuracy"],
                "importance": "5-7% JEE weightage. Easy marks — always appears in JEE Main. Master dimensional analysis for formula derivation."
            },
            "Vectors": {
                "difficulty": "Medium",
                "description": "Vectors are the mathematical backbone of all of physics. Unlike scalars (like temperature or mass), vectors carry both magnitude and direction — making them essential for describing forces, velocities, electric fields, and more. This chapter covers vector addition, subtraction, dot product, cross product, and resolution of vectors into components. Mastery of vectors is non-negotiable for JEE since they appear in mechanics, electromagnetism, and optics.",
                "roadmap": ["Scalars vs Vectors", "Vector Addition (Triangle & Parallelogram Law)", "Resolution of Vectors", "Dot Product (Scalar Product)", "Cross Product (Vector Product)"],
                "topics": ["Unit Vectors î, ĵ, k̂", "Vector Addition & Subtraction", "Component Method of Addition", "Scalar Product & its properties", "Vector Product & Right-Hand Rule", "Position & Displacement Vectors", "Angle between two vectors"],
                "formulas": ["A·B = AB cosθ (Dot Product)", "|A×B| = AB sinθ (Cross Product)", "R² = A² + B² + 2AB cosθ (Resultant)", "tan α = B sinθ/(A + B cosθ) (Direction)", "Magnitude = √(Ax² + Ay² + Az²)"],
                "laws": ["Triangle Law of Vector Addition", "Parallelogram Law of Vector Addition"],
                "applications": ["Finding resultant force, velocity", "Torque calculation (r × F)", "Work done (F · displacement)"],
                "common_mistakes": ["Confusing dot and cross product results", "Not using right-hand rule for cross product", "Adding magnitudes instead of vectors"],
                "importance": "Used in every chapter of physics. 2-3 direct JEE questions. Essential tool for the entire syllabus."
            },
            "Motion in One Dimension": {
                "difficulty": "Easy",
                "description": "Motion in One Dimension (Kinematics in 1D) is the study of objects moving along a straight line. It introduces the core concepts of displacement, velocity, and acceleration with their scalar and vector nature understood in a single dimension. The three equations of motion derived from uniform acceleration are the most-used formulas in all of physics. Velocity-time and displacement-time graphs are powerful visual tools for analysing motion and are frequently tested in JEE Main.",
                "roadmap": ["Distance vs Displacement", "Speed vs Velocity", "Uniform & Non-uniform Acceleration", "Equations of Motion", "Graphs: s-t, v-t, a-t"],
                "topics": ["Average & Instantaneous Velocity", "Uniform Acceleration Equations", "Free Fall & Motion Under Gravity", "Relative Motion in 1D", "Analysis of v-t graphs for displacement", "Reaction Time"],
                "formulas": ["v = u + at", "s = ut + ½at²", "v² = u² + 2as", "s_nth = u + a(2n−1)/2", "g = 9.8 m/s² (free fall)"],
                "applications": ["Free fall problems", "Particle motion analysis from graphs", "Relative velocity between trains/cars"],
                "common_mistakes": ["Not taking sign convention consistently", "Confusing distance and displacement in free fall", "Misreading velocity-time graphs"],
                "importance": "Foundation of mechanics. 1-2 questions in JEE Main. Easy scoring with graph-based problems."
            },
            "Motion in Two Dimensions": {
                "difficulty": "Medium",
                "description": "Motion in Two Dimensions extends kinematics to a plane, where an object moves under the influence of forces in both x and y directions independently. The most important application is Projectile Motion — the path of a ball thrown at an angle, which traces a parabola. Uniform Circular Motion introduces centripetal acceleration and is the base for topics like banking of roads, conical pendulums, and planetary orbits. Relative velocity in 2D (river-boat problems, rain-man problems) is a JEE favourite.",
                "roadmap": ["2D Kinematics — Component Method", "Projectile Motion (Horizontal & Angled)", "Horizontal Projectile", "Uniform Circular Motion", "Relative Velocity in 2D"],
                "topics": ["Projectile Range, Height, Time of Flight", "Equation of Trajectory y = x tanθ − gx²/2u²cos²θ", "Centripetal & Centrifugal Acceleration", "Banking of Roads", "River-Boat & Rain-Man Problems", "Conical Pendulum"],
                "formulas": ["R = u²sin2θ/g", "H = u²sin²θ/2g", "T = 2u sinθ/g", "aₖ = v²/r = ω²r (centripetal)", "tan θ = v²/rg (banking)"],
                "applications": ["Ball throw, cannon fire (projectile)", "Satellite circular motion", "Car on banked curves"],
                "common_mistakes": ["Not resolving initial velocity into components", "Using g = 10 inconsistently", "Forgetting centripetal is always towards center"],
                "importance": "High weightage — 2-3 JEE questions. Projectile motion is one of the most tested topics in JEE Main."
            },
            "Laws of Motion": {
                "difficulty": "Hard",
                "description": "Newton's Laws of Motion are the cornerstone of classical mechanics. The First Law defines inertia — a body at rest or uniform motion stays so unless acted upon by a net external force. The Second Law (F = ma) quantifies force. The Third Law states every action has an equal and opposite reaction. This chapter involves complex systems: pulley problems, connected bodies, friction on inclined planes, circular motion dynamics, and pseudo forces in non-inertial frames — all JEE favourites.",
                "roadmap": ["Newton's Three Laws", "Free Body Diagram (FBD) Technique", "Friction — Static, Kinetic, Rolling", "Connected Bodies & Pulleys", "Pseudo Force in Non-Inertial Frames"],
                "topics": ["Inertia & Mass", "Normal Force", "Tension in Strings", "Static & Kinetic Friction", "Angle of Friction & Angle of Repose", "Banking Without Friction", "Atwood Machine", "Wedge Problems"],
                "formulas": ["F = ma", "p = mv", "Impulse J = FΔt = Δp", "f_s ≤ μ_s N", "f_k = μ_k N", "tan φ = μ (angle of friction)"],
                "laws": ["Newton's 1st Law — Law of Inertia", "Newton's 2nd Law — F = dp/dt", "Newton's 3rd Law — Action-Reaction"],
                "applications": ["Car braking distance (friction)", "Lift problems (apparent weight)", "Rocket propulsion (3rd law)"],
                "common_mistakes": ["Drawing incomplete free body diagrams", "Applying pseudo force in inertial frame", "Taking friction direction wrong"],
                "importance": "15-20% of JEE mechanics questions. Mastering FBD is the single most important skill for JEE physics."
            },
            "Work, Energy & Power": {
                "difficulty": "Medium",
                "description": "Work, Energy & Power is one of the most elegant chapters in physics. Work is done when a force produces displacement in its direction. Kinetic energy is the energy of motion; potential energy is stored energy. The Work-Energy Theorem directly connects net work to change in kinetic energy. Conservation of mechanical energy is a powerful principle that bypasses force analysis entirely. Elastic and inelastic collisions and the concept of power round off this JEE-critical chapter.",
                "roadmap": ["Work by Constant & Variable Force", "Kinetic & Potential Energy", "Work-Energy Theorem", "Conservation of Mechanical Energy", "Power & Efficiency"],
                "topics": ["Work done at an angle (W = Fs cosθ)", "Work done by spring force", "Conservative & Non-Conservative Forces", "Potential Energy Curve", "Elastic & Inelastic Collisions", "Coefficient of Restitution", "Instantaneous Power"],
                "formulas": ["W = F·s·cosθ", "W = ½kx² (spring)", "KE = ½mv²", "PE = mgh", "W_net = ΔKE", "P = F·v", "e = relative speed of separation/relative speed of approach"],
                "applications": ["Roller coaster (energy conservation)", "Spring systems", "Car engine power output"],
                "common_mistakes": ["Using W = Fs without cosθ", "Not accounting for friction work in energy conservation", "Confusing elastic collision (KE conserved) with inelastic"],
                "importance": "10-12% JEE weightage. Conservation of energy is the most powerful problem-solving tool in mechanics."
            },
            "Centre of Mass, Momentum & Collisions": {
                "difficulty": "Hard",
                "description": "The Centre of Mass (COM) is the single point that represents the entire mass distribution of a system. For external forces, the entire system behaves as if all mass is concentrated at the COM. Conservation of linear momentum — one of the deepest laws in physics — holds when no external force acts. Elastic collisions conserve both momentum and kinetic energy; inelastic do not. This chapter is essential for solving complex JEE problems involving explosions, rocket propulsion, and oblique collisions.",
                "roadmap": ["Centre of Mass — Discrete & Continuous Systems", "Motion of COM", "Conservation of Linear Momentum", "Types of Collisions", "Coefficient of Restitution"],
                "topics": ["COM of triangles, semicircles, hemispheres", "Velocity & Acceleration of COM", "Explosion problems", "Perfectly Elastic Collision formulas", "Perfectly Inelastic Collision", "Oblique Collisions", "Rocket Propulsion (Tsiolkovsky equation)"],
                "formulas": ["x_COM = Σmᵢxᵢ / Σmᵢ", "p_total = constant (no ext. force)", "In elastic: v₁' = (m₁-m₂)u₁/(m₁+m₂) + 2m₂u₂/(m₁+m₂)", "e = 1 (elastic), e = 0 (perfectly inelastic)", "v_rocket = u ln(M₀/M)"],
                "applications": ["Explosion breaking a body apart", "Ballistic pendulum", "Rocket thrust"],
                "common_mistakes": ["Applying momentum conservation when external forces exist", "Forgetting COM velocity during explosion = 0 change"],
                "importance": "8-10% JEE weightage. Collisions appear in almost every JEE paper."
            },
            "Rotational Motion": {
                "difficulty": "Hard",
                "description": "Rotational Motion is the most mathematically rich chapter in Class 11 Physics. Every translational concept has a rotational analogue: force → torque, mass → moment of inertia, linear momentum → angular momentum. The parallel and perpendicular axis theorems allow calculation of moment of inertia for complex bodies. Rolling motion combines both translational and rotational KE. Angular momentum conservation (like a figure skater pulling arms in) is one of the most beautiful applications in physics.",
                "roadmap": ["Angular Kinematics (α, ω, θ)", "Torque & Couple", "Moment of Inertia & Theorems", "Angular Momentum & Conservation", "Rolling Motion Without Slipping"],
                "topics": ["Equations of Rotational Motion", "Torque = r × F", "MI of rod, disk, ring, sphere, hollow sphere", "Parallel Axis Theorem: I = I_cm + Md²", "Perpendicular Axis Theorem: I_z = I_x + I_y", "Rolling without slipping condition v = Rω", "Angular Impulse"],
                "formulas": ["τ = Iα", "L = Iω", "KE_roll = ½Iω² + ½mv²", "I_disk = ½MR²", "I_sphere = ⅖MR²", "I_ring = MR²", "I_rod(center) = ML²/12"],
                "laws": ["Conservation of Angular Momentum: L = constant when τ_ext = 0"],
                "applications": ["Spinning top & gyroscope", "Figure skater pulling arms in", "Rolling objects down inclines"],
                "common_mistakes": ["Using Iα without checking rolling constraint", "Forgetting to add ½mv² for rolling KE", "Wrong axis for parallel axis theorem"],
                "importance": "15% JEE weightage. One of the highest-scoring chapters in JEE Advanced. Master MI of standard bodies."
            },
            "Gravitation": {
                "difficulty":": Medium",
                "description": "Gravitation is the fundamental force that governs the motion of planets, moons, and satellites. Newton's Law of Universal Gravitation gives the attractive force between any two masses. This chapter covers gravitational field, potential, escape velocity, orbital velocity, geostationary satellites, and Kepler's three laws. Variations of 'g' with altitude, depth, latitude, and Earth's rotation are tested frequently in JEE Main. Understanding gravitational potential energy well is key to solving orbital mechanics problems.",
                "roadmap": ["Newton's Law of Gravitation", "Gravitational Field & Potential", "Variation of g", "Orbital & Escape Velocity", "Satellites & Kepler's Laws"],
                "topics": ["Gravitational Field Intensity g = F/m", "Gravitational Potential V = −GM/r", "Variation of g with altitude: g' = g(1−2h/R)", "Variation of g with depth: g' = g(1−d/R)", "Orbital velocity v₀ = √(GM/r)", "Escape velocity vₑ = √(2GM/R)", "Geostationary Satellites (T = 24 hr)", "Kepler's 3 Laws"],
                "formulas": ["F = Gm₁m₂/r²", "g = GM/R²", "vₑ = √(2gR) = 11.2 km/s", "v₀ = √(gR) = 7.9 km/s", "T² ∝ r³ (Kepler's 3rd)", "PE = −GMm/r"],
                "laws": ["Kepler's 1st: Elliptical orbits", "Kepler's 2nd: Equal areas in equal times", "Kepler's 3rd: T² ∝ a³"],
                "applications": ["Moon's orbital period", "GPS & communication satellites", "Space mission trajectory planning"],
                "common_mistakes": ["Confusing g at surface vs. at height h", "Using vₑ = √(2gR) only at Earth's surface"],
                "importance": "8-10% JEE weightage. Numericals on satellites, escape velocity appear every year."
            },
            "Mechanical Properties of Solids": {
                "difficulty": "Medium",
                "description": "Mechanical Properties of Solids explores how solid materials respond to external forces — whether they stretch, compress, or shear. Stress is the internal restoring force per unit area; strain is the fractional deformation. The three elastic moduli (Young's, Bulk, Shear) quantify this response. Hooke's Law states that stress is proportional to strain within the elastic limit. This chapter has direct applications in engineering — from bridges to aircraft wings — and appears in JEE in both conceptual and numerical forms.",
                "roadmap": ["Stress & Strain", "Hooke's Law & Elastic Limit", "Young's Modulus", "Bulk & Shear Modulus", "Elastic Potential Energy"],
                "topics": ["Tensile, Compressive & Shear Stress", "Longitudinal, Volumetric & Shear Strain", "Young's Modulus Y = Stress/Strain", "Bulk Modulus B = −P/(ΔV/V)", "Modulus of Rigidity η = Shear stress/Shear strain", "Poisson's Ratio σ = −(lateral strain)/(longitudinal strain)", "Elastic PE = ½ × Stress × Strain × Volume"],
                "formulas": ["Y = FL/AΔL", "B = −VΔP/ΔV", "Elastic PE = ½YAε²l", "Work done in stretching = ½ × F × ΔL"],
                "applications": ["Steel vs rubber elasticity", "Bridges and building materials", "Bone mechanics"],
                "common_mistakes": ["Confusing modulus of rigidity with Young's modulus", "Using wrong formula for elastic PE"],
                "importance": "3-5% JEE weightage. Short but conceptual. Know all 3 moduli and Hooke's law region."
            },
            "Mechanical Properties of Fluids": {
                "difficulty": "Medium",
                "description": "Mechanical Properties of Fluids bridges the gap between solid mechanics and fluid dynamics. Pressure in a fluid acts in all directions and increases with depth. Archimedes' Principle explains buoyancy — the basis of ships and submarines. Bernoulli's Theorem (conservation of energy in fluid flow) explains aircraft lift, venturimeter, and why a shower curtain blows inward. Viscosity and surface tension are real-world phenomena governing everyday experiences from blood flow to soap bubbles.",
                "roadmap": ["Fluid Pressure & Pascal's Law", "Archimedes' Principle & Buoyancy", "Bernoulli's Theorem", "Viscosity & Stoke's Law", "Surface Tension & Capillarity"],
                "topics": ["Hydrostatic Pressure P = ρgh", "Pascal's Law & Hydraulic Machines", "Buoyant Force = ρ_fluid × V_submerged × g", "Equation of Continuity A₁v₁ = A₂v₂", "Torricelli's theorem", "Venturimeter & Pitot Tube", "Stoke's Law F = 6πηrv", "Terminal Velocity", "Excess Pressure in bubble/drop", "Capillary Rise h = 2T cosθ/ρrg"],
                "formulas": ["Bernoulli: P + ½ρv² + ρgh = constant", "h = 2T cosθ/(ρrg)", "v_terminal = 2r²(ρ-σ)g/9η", "Excess pressure inside bubble = 4T/r", "Inside drop = 2T/r"],
                "applications": ["Aircraft wing lift (Bernoulli)", "Hydraulic brakes/lift", "Blood flow through arteries"],
                "common_mistakes": ["Applying Bernoulli where flow is not streamline", "Forgetting factor of 2T vs 4T for drops vs bubbles"],
                "importance": "8-10% JEE weightage. Bernoulli and surface tension are frequently tested with tricky numericals."
            },
            "Thermal Properties of Matter & Heat Transfer": {
                "difficulty": "Medium",
                "description": "This chapter covers the fundamentals of heat, temperature, and how heat flows between objects. Temperature scales (Celsius, Fahrenheit, Kelvin) and their conversions are basic. Thermal expansion — linear, area, and volumetric — explains why railway tracks have gaps and thermometers work. Specific heat capacity and latent heat govern heating and phase changes. The three modes of heat transfer — conduction (metals), convection (fluids), and radiation (electromagnetic waves) — are all tested in JEE with both conceptual and calculation-based questions.",
                "roadmap": ["Temperature Scales & Conversion", "Thermal Expansion", "Specific Heat & Calorimetry", "Latent Heat & Phase Changes", "Modes of Heat Transfer"],
                "topics": ["Linear: ΔL = αLΔT, Area: ΔA = 2αAΔT, Volume: ΔV = 3αVΔT", "Specific Heat Q = mcΔT", "Latent Heat Q = mL", "Calorimetry — Heat gained = Heat lost", "Newton's Law of Cooling", "Stefan-Boltzmann Law E = σT⁴", "Wien's Displacement Law λ_max T = b", "Thermal Resistance & Conduction Rate"],
                "formulas": ["Q = mcΔT", "Q = mL", "dQ/dt = kA(ΔT/x) (conduction)", "E = εσT⁴ (radiation)", "λ_max = b/T = 0.29 cm·K/T"],
                "laws": ["Newton's Law of Cooling", "Stefan-Boltzmann Law", "Wien's Displacement Law"],
                "applications": ["Thermos flask (radiation prevention)", "Land & sea breezes (convection)", "Pressure cooker (latent heat)"],
                "common_mistakes": ["Confusing specific heat with heat capacity", "Not converting temperature to Kelvin for radiation formulas"],
                "importance": "8-10% JEE weightage. Wien's law and Stefan's law appear regularly in JEE Main."
            },
            "Thermodynamics": {
                "difficulty": "Hard",
                "description": "Thermodynamics is the science of energy transformation. It studies how heat and work are interconverted and what limits exist on such conversions. The Zeroth Law defines temperature; the First Law is conservation of energy for heat engines; the Second Law introduces entropy and the concept of irreversibility. The Carnot engine represents the theoretical maximum efficiency of any heat engine. P-V diagrams (isothermal, adiabatic, isochoric, isobaric processes) are core JEE exam tools.",
                "roadmap": ["Thermodynamic Systems & Processes", "Zeroth & First Laws", "PV Diagrams — Isothermal, Adiabatic, Isochoric, Isobaric", "Second Law & Carnot Engine", "Entropy"],
                "topics": ["Internal Energy ΔU", "Work done by gas W = ∫PdV", "Isothermal: T constant, W = nRT ln(V₂/V₁)", "Adiabatic: PVγ = constant", "Isochoric: W = 0", "Isobaric: W = PΔV", "Carnot Cycle — 4 steps", "Efficiency η = 1 − T_cold/T_hot"],
                "formulas": ["ΔU = Q − W (First Law)", "W = PΔV (isobaric)", "W = nRT ln(V₂/V₁) (isothermal)", "PVγ = const (adiabatic)", "η_Carnot = 1 − T₂/T₁", "COP (refrigerator) = T₂/(T₁−T₂)"],
                "laws": ["Zeroth Law — Thermal Equilibrium defines temperature", "First Law — Energy Conservation", "Second Law — Heat flows from hot to cold spontaneously"],
                "applications": ["Car engines (heat engine)", "Refrigerators & ACs (reverse heat engine)", "Steam turbines"],
                "common_mistakes": ["Forgetting sign convention: Q positive when absorbed, W positive when done by gas", "Using adiabatic formula for isothermal process"],
                "importance": "10-12% JEE weightage. Carnot efficiency and PV diagrams appear in almost every JEE paper."
            },
            "Kinetic Theory of Gases": {
                "difficulty": "Medium",
                "description": "The Kinetic Theory of Gases provides the microscopic explanation for the macroscopic behaviour described by gas laws. It models gas as a large number of tiny particles in random motion, colliding elastically. From this model, we can derive pressure, temperature (a measure of average KE), and the ideal gas equation PV = nRT. The Maxwell speed distribution gives us RMS, average, and most probable speeds. Degrees of freedom and the equipartition theorem explain the heat capacity of gases.",
                "roadmap": ["Postulates of Kinetic Theory", "Pressure & Temperature from Molecular Model", "Gas Laws — Boyle's, Charles's, Gay-Lussac's", "RMS, Average & Most Probable Speed", "Degrees of Freedom & Equipartition Theorem"],
                "topics": ["Ideal Gas Equation PV = nRT", "KE of gas molecule = 3/2 kT", "v_rms = √(3RT/M)", "v_avg = √(8RT/πM)", "v_mp = √(2RT/M)", "Mean Free Path λ = 1/(√2 nπd²)", "Degrees of freedom: monoatomic=3, diatomic=5, polyatomic=6", "Cv = fR/2, Cp = (f+2)R/2, γ = Cp/Cv"],
                "formulas": ["PV = nRT", "KE = ½mv² = 3/2 kT per molecule", "v_rms : v_avg : v_mp = √3 : √(8/π) : √2", "γ = 5/3 (monoatomic), 7/5 (diatomic)"],
                "applications": ["Explaining Boyle's and Charles's laws from molecular model", "Gas thermometers", "Effusion and diffusion of gases"],
                "common_mistakes": ["Confusing v_rms, v_avg, and v_mp", "Using wrong degrees of freedom for diatomic vs monoatomic"],
                "importance": "8% JEE weightage. Speed ratios and Cv, Cp, γ values are must-know formulas."
            },
            "Waves & Sound": {
                "difficulty": "Medium",
                "description": "Waves chapter explores how energy propagates through a medium as a disturbance. Waves are classified as transverse (light) or longitudinal (sound). The wave equation y = A sin(kx − ωt) encodes amplitude, frequency, and speed. Superposition of waves creates interference (constructive/destructive), beats, and standing waves (stationary waves). Resonance in strings and air columns explains musical instruments. Doppler effect explains the change in frequency with relative motion — from ambulance sirens to red-shift of galaxies.",
                "roadmap": ["Wave Parameters — Amplitude, Frequency, Wavelength, Speed", "Wave Equation & Phase", "Superposition & Interference", "Standing Waves — Strings & Pipes", "Beats & Doppler Effect"],
                "topics": ["Transverse vs Longitudinal Waves", "Wave Speed v = fλ", "Phase difference vs Path difference: Δφ = (2π/λ)Δx", "Nodes & Antinodes in standing waves", "Harmonics & Overtones", "Open & Closed Pipe frequencies", "Beat frequency = |f₁ − f₂|", "Doppler: f' = f(v±v_observer)/(v∓v_source)"],
                "formulas": ["v = fλ", "y = A sin(kx − ωt)", "v_sound = √(γP/ρ) = √(γRT/M)", "f_n (string) = n/2L √(T/μ)", "f_n (open pipe) = nv/2L", "f_n (closed pipe) = (2n−1)v/4L"],
                "laws": ["Principle of Superposition"],
                "applications": ["Music (string instruments, wind instruments)", "Medical ultrasound imaging", "Speed gun (Doppler radar)"],
                "common_mistakes": ["Forgetting closed pipe has only odd harmonics", "Wrong sign in Doppler formula for approaching vs receding"],
                "importance": "8-10% JEE weightage. Standing waves, Doppler effect are annual JEE favourites. Must practise both strings and pipes."
            },
            "Simple Harmonic Motion": {
                "difficulty": "Medium",
                "description": "Simple Harmonic Motion (SHM) is the most important kind of oscillatory motion in physics. It occurs when a restoring force acts proportional to displacement: F = −kx. The motion is sinusoidal in time, with well-defined amplitude, time period, frequency, and phase. Spring-mass systems and simple pendulums are the two key oscillators. Energy in SHM continuously exchanges between KE and PE, with total energy conserved. SHM is the foundation for understanding LC circuits, resonance, and wave motion.",
                "roadmap": ["Conditions for SHM", "Displacement, Velocity, Acceleration in SHM", "Spring-Mass System", "Simple & Compound Pendulum", "Energy in SHM"],
                "topics": ["F = −kx (restoring force)", "x = A sin(ωt + φ)", "v = Aω cos(ωt + φ)", "a = −ω²x", "T = 2π√(m/k) — spring", "T = 2π√(l/g) — pendulum", "KE = ½mω²(A²−x²)", "PE = ½mω²x²", "Total E = ½mω²A²", "Damped & Forced Oscillations"],
                "formulas": ["x = A sin(ωt + φ)", "T = 2π/ω", "T_spring = 2π√(m/k)", "T_pendulum = 2π√(L/g)", "E_total = ½kA²"],
                "applications": ["Pendulum clock", "Car suspension (damped oscillations)", "Resonance bridges (Tacoma Narrows)"],
                "common_mistakes": ["Confusing time period with frequency", "Forgetting amplitude A is max displacement from equilibrium", "Using T = 2π√(l/g) for all pendulums — valid only for small angles"],
                "importance": "10-12% JEE weightage. SHM is combined with waves, circular motion, and LC circuits in advanced problems."
            }
        },

        "Chemistry": {
            "Some Basic Concepts of Chemistry": {
                "difficulty": "Easy",
                "description": "Some Basic Concepts of Chemistry lays the quantitative foundation of chemistry. It covers Dalton's Atomic Theory, laws of chemical combination (conservation of mass, definite proportions, multiple proportions), SI units, significant figures, mole concept, molar mass, and stoichiometry. Every chemistry calculation — from titrations to gas laws — builds on this chapter. The mole is the single most important concept: it bridges the atomic world (atoms, molecules) with the macroscopic world (grams, litres).",
                "roadmap": ["Matter & its Nature", "Laws of Chemical Combination", "Atomic & Molecular Masses", "Mole Concept & Molar Mass", "Stoichiometry & Limiting Reagent"],
                "topics": ["Dalton's Atomic Theory", "Law of Conservation of Mass", "Law of Definite Proportions", "Law of Multiple Proportions", "Avogadro's Law", "Atomic mass unit (amu)", "Mole = 6.022×10²³ particles", "Molar mass", "Percentage composition", "Empirical & Molecular formula", "Limiting reagent", "% yield"],
                "formulas": ["Moles = mass/molar mass", "Moles = volume(L)/22.4 (at STP)", "% composition = (mass of element / molar mass) × 100", "Empirical formula from % composition"],
                "applications": ["Titration calculations", "Finding limiting reagent", "Industrial chemical production"],
                "common_mistakes": ["Confusing empirical and molecular formula", "Forgetting to convert g to mol before stoichiometry", "Wrong significant figures"],
                "importance": "8-10% JEE Main weightage. Easy marks — mole concept and stoichiometry appear every year."
            },
            "Structure of Atom": {
                "difficulty": "Hard",
                "description": "Structure of Atom traces the development of our understanding of atomic architecture. Thomson's plum pudding model gave way to Rutherford's nuclear model, then Bohr quantised the hydrogen atom. The quantum mechanical model — built on de Broglie's wave-particle duality and Heisenberg's uncertainty principle — provides the modern picture. Quantum numbers (n, l, ml, ms) define each electron's state, and orbital shapes (s, p, d) determine chemical bonding.",
                "roadmap": ["Thomson & Rutherford Models", "Bohr's Model of Hydrogen", "Dual Nature & Uncertainty Principle", "Quantum Mechanical Model", "Electronic Configuration"],
                "topics": ["Cathode rays & discovery of electron", "Rutherford's gold foil experiment", "Atomic number Z, mass number A", "Bohr's postulates: rn = 0.529 n² Å, En = -13.6/n² eV", "Hydrogen spectrum — Lyman, Balmer, Paschen series", "de Broglie: λ = h/mv", "Heisenberg: Δx·Δp ≥ h/4π", "Quantum numbers: n, l, ml, ms", "Shapes of s, p, d orbitals", "Aufbau principle, Pauli's exclusion, Hund's rule", "Extra stability of half-filled and fully filled orbitals"],
                "formulas": ["E = hν", "λ = h/p = h/mv (de Broglie)", "Δx·Δp ≥ h/4π (Heisenberg)", "1/λ = R(1/n₁² - 1/n₂²) (Rydberg)", "En = -13.6/n² eV (Bohr, H atom)"],
                "applications": ["Flame tests (atomic spectra)", "Laser technology", "Electron microscopes"],
                "common_mistakes": ["Forgetting max electrons per orbital = 2", "Wrong orbital filling order (remember: 4s before 3d)", "Confusing quantum numbers l and ml"],
                "importance": "10-12% JEE weightage. Electronic config and quantum numbers appear in every JEE paper."
            },
            "Classification of Elements & Periodicity": {
                "difficulty": "Medium",
                "description": "Classification of Elements & Periodicity explains the logic behind the periodic table — how Mendeleev arranged 63 known elements by atomic mass, and how Moseley's discovery of atomic number led to the modern periodic law. The chapter covers periodic trends: atomic radius, ionic radius, ionization energy, electron affinity, electronegativity, and metallic character. These trends arise from nuclear charge and electron shielding and predict reactivity patterns across groups and periods.",
                "roadmap": ["Historical Development of Periodic Table", "Modern Periodic Law", "Blocks (s, p, d, f)", "Periodic Trends — Atomic Radius", "Periodic Trends — IE, EA, EN"],
                "topics": ["Mendeleev's periodic law", "Modern periodic law: properties are periodic functions of atomic number", "s, p, d, f block elements", "Atomic radius trends: increases down group, decreases across period", "Ionic radius: cation < parent atom < anion", "Ionization energy: IE₁ < IE₂ < IE₃", "Electron affinity", "Electronegativity: Pauling scale", "Metallic & non-metallic character", "Diagonal relationship (Li-Mg, Be-Al)"],
                "formulas": ["Effective nuclear charge (Zeff) = Z - σ (Slater's rules)", "IE trend across period increases", "Atomic radius: increases down group"],
                "applications": ["Predicting compound formation", "Explaining anomalous first ionization energies", "Diagonal relationships in industry"],
                "common_mistakes": ["Forgetting EA is exothermic (negative)", "Wrong trend for electron affinity (not always regular)", "Confusing atomic and ionic radius trends"],
                "importance": "8-10% JEE weightage. Mostly conceptual — easy marks with good memory."
            },
            "Chemical Bonding & Molecular Structure": {
                "difficulty": "Hard",
                "description": "Chemical Bonding & Molecular Structure explains why and how atoms bond. Ionic bonds form by electron transfer; covalent bonds by sharing. VSEPR theory predicts molecular geometry from electron pair repulsion. Valence Bond Theory explains hybridization (sp, sp², sp³, sp³d, sp³d²). Molecular Orbital Theory explains paramagnetism and bond order in diatomic molecules. Fajan's rules, hydrogen bonding, and dipole moment are essential IIT JEE topics.",
                "roadmap": ["Octet Rule & Exceptions", "Ionic vs Covalent Bonding", "Lewis Dot Structures", "VSEPR Theory — Molecular Geometry", "Hybridization & MO Theory"],
                "topics": ["Electronegativity & Fajan's rules", "Lewis dot structures", "Formal charge: FC = V - L - B/2", "VSEPR — shapes: linear, bent, trigonal planar, tetrahedral, trigonal bipyramidal, octahedral", "Hybridization: sp (linear), sp² (planar), sp³ (tetrahedral)", "sp³d (trigonal bipyramidal), sp³d² (octahedral)", "Bond angle variations (lone pairs reduce angles)", "Pi and sigma bonds", "MO Theory: bonding and antibonding orbitals", "Bond order = (Nb - Na)/2", "Paramagnetism of O₂ (explained by MO)", "Hydrogen bonding — intra and inter molecular", "Dipole moment"],
                "formulas": ["Bond order = (Nb - Na)/2", "FC = V - L - B/2", "% Ionic character = 16|Δχ| + 3.5(Δχ)²"],
                "applications": ["Predicting molecular shapes", "Explaining boiling points via H-bonding", "DNA double helix (H-bonds)"],
                "common_mistakes": ["Wrong number of lone pairs in VSEPR", "Confusing sigma/pi with single/double bonds", "Forgetting O₂ is paramagnetic"],
                "importance": "12-15% JEE weightage. VSEPR and MO theory problems appear every year."
            },
            "States of Matter": {
                "difficulty": "Medium",
                "description": "States of Matter explores how matter behaves as a gas, liquid, or solid under different conditions of temperature and pressure. The gaseous state is governed by gas laws — Boyle's, Charles's, Graham's, Dalton's, and the Ideal Gas Law PV = nRT. Kinetic molecular theory explains microscopic reasons behind macroscopic properties. Real gas deviations from ideality are explained by van der Waals equation. The liquid state introduces surface tension, viscosity, and vapour pressure.",
                "roadmap": ["Gas Laws — Boyle, Charles, Graham", "Ideal Gas Equation", "Dalton's Law of Partial Pressure", "Kinetic Theory of Gases", "Real Gases & van der Waals"],
                "topics": ["Boyle's Law: PV = constant (constant T)", "Charles's Law: V/T = constant (constant P)", "Avogadro's Law: V ∝ n (constant T, P)", "Graham's Law: r₁/r₂ = √(M₂/M₁)", "Dalton's Law: P_total = P₁ + P₂ + P₃...", "PV = nRT (Ideal Gas)", "Kinetic energy = 3/2 kT per molecule", "Compressibility factor Z = PV/nRT", "van der Waals equation: (P + an²/V²)(V - nb) = nRT", "Liquefaction of gases — Joule-Thomson effect", "Surface tension, viscosity, vapour pressure"],
                "formulas": ["PV = nRT", "r₁/r₂ = √(M₂/M₁)", "(P + an²/V²)(V - nb) = nRT", "Z = PV/nRT", "vrms = √(3RT/M)"],
                "applications": ["Weather balloons", "Scuba diving (Dalton's Law)", "Industrial gas storage"],
                "common_mistakes": ["Using wrong R value (8.314 vs 0.0821 L·atm/mol·K)", "Forgetting to convert temperature to Kelvin"],
                "importance": "8-10% JEE weightage. Gas law numericals are straightforward — easy marks."
            },
            "Thermodynamics (Chemistry)": {
                "difficulty": "Hard",
                "description": "Chemical Thermodynamics studies energy changes in chemical processes. The First Law (energy conservation) introduces internal energy (U), heat (q), and work (w). Hess's Law allows calculation of reaction enthalpies indirectly. Entropy (S) measures disorder, and the Second Law defines spontaneity. Gibbs Free Energy (G = H - TS) determines whether a reaction will occur spontaneously — a central concept for chemical equilibrium and electrochemistry.",
                "roadmap": ["System, Surroundings & State Functions", "First Law — Internal Energy & Enthalpy", "Hess's Law", "Second Law — Entropy", "Gibbs Free Energy & Spontaneity"],
                "topics": ["System, surroundings, boundary", "Open, closed, isolated systems", "State functions: U, H, S, G", "First Law: ΔU = q + w", "Work: w = -PΔV (expansion)", "Enthalpy: ΔH = ΔU + PΔV", "Exothermic (ΔH < 0) & Endothermic (ΔH > 0)", "Standard enthalpy of formation, combustion, neutralization, bond dissociation", "Hess's Law: ΔHrxn = ΣΔHproducts - ΣΔHreactants", "Entropy S: ΔS > 0 for spontaneous at constant U", "Gibbs: ΔG = ΔH - TΔS", "ΔG° = -RT ln K = -nFE°"],
                "formulas": ["ΔU = q + w", "ΔH = ΔU + ΔnRT (for gases)", "ΔG = ΔH - TΔS", "ΔG° = -RT ln K", "ΔG° = -nFE°"],
                "applications": ["Predicting reaction spontaneity", "Fuel cell efficiency", "Industrial process design"],
                "common_mistakes": ["Sign error: work done ON system is +w", "Forgetting ΔG < 0 for spontaneous reaction", "Wrong sign for enthalpy of combustion"],
                "importance": "12-15% JEE weightage. High scoring in JEE Advanced — ΔG = -RT ln K is critical."
            },
            "Equilibrium": {
                "difficulty": "Hard",
                "description": "Equilibrium is one of the most concept-heavy chapters in JEE Chemistry. Chemical equilibrium explains why reactions don't go to completion — forward and reverse reactions occur simultaneously. Kc and Kp quantify the equilibrium position. Le Chatelier's Principle predicts how equilibrium shifts with changes in concentration, pressure, temperature, and catalyst. Ionic equilibrium covers weak acid/base dissociation, pH, buffer solutions, solubility product Ksp, and common ion effect.",
                "roadmap": ["Law of Chemical Equilibrium", "Kc & Kp", "Le Chatelier's Principle", "Ionic Equilibrium — Acids & Bases", "Buffer & Solubility"],
                "topics": ["Dynamic equilibrium", "Equilibrium constant expression", "Kp = Kc(RT)^Δn", "Reaction quotient Q: Q < K (forward), Q > K (reverse)", "Le Chatelier: effect of concentration, pressure, temperature, catalyst", "Brønsted-Lowry acids & bases", "Lewis acids & bases", "Conjugate acid-base pairs", "Ka, Kb, Kw = Ka × Kb = 10⁻¹⁴", "pH = -log[H⁺]", "Buffer: pH = pKa + log([A⁻]/[HA])", "Ksp and solubility", "Common ion effect"],
                "formulas": ["Kp = Kc(RT)^Δn", "pH = -log[H⁺]", "Kw = [H⁺][OH⁻] = 10⁻¹⁴", "Henderson-Hasselbalch: pH = pKa + log([A⁻]/[HA])", "Ksp = [A⁺][B⁻] for AB"],
                "applications": ["Industrial Haber process (NH₃ synthesis)", "Buffer in blood (carbonate buffer pH 7.4)", "Stalagmites and stalactites (Ksp)"],
                "common_mistakes": ["Including solids/liquids in Kc/Kp expressions", "Confusing Q with K", "Wrong buffer pH formula"],
                "importance": "15-18% JEE weightage. Maximum JEE questions come from this chapter."
            },
            "Redox Reactions": {
                "difficulty": "Medium",
                "description": "Redox Reactions (Reduction-Oxidation) involve electron transfer between chemical species. Oxidation is loss of electrons (OIL — Oxidation Is Loss), reduction is gain (GIG — Gain Is Gain). The concept of oxidation numbers allows identification of redox changes and balancing complex equations using the half-reaction method or oxidation number method. Understanding redox is essential for electrochemistry, metallurgy, and biological processes.",
                "roadmap": ["Oxidation & Reduction Definitions", "Oxidation Number Rules", "Identifying Redox Reactions", "Balancing Redox Equations — Half-Reaction Method", "Disproportionation Reactions"],
                "topics": ["OIL-RIG mnemonic", "Oxidation number rules (O = -2, H = +1 etc.)", "Half-reaction method in acidic/basic medium", "Disproportionation: same element oxidized & reduced", "Common oxidizing agents: KMnO₄, K₂Cr₂O₇, H₂O₂, HNO₃", "Common reducing agents: Fe, C, Na₂S₂O₃"],
                "formulas": ["Electrons lost = electrons gained", "Net ionic equation method"],
                "applications": ["Corrosion of metals (rusting)", "Bleaching powder action", "Metabolic reactions (respiration)"],
                "common_mistakes": ["Wrong oxidation number for S, N, Cl (variable valency)", "Forgetting to balance charge with H⁺ or OH⁻", "Misidentifying oxidizing vs reducing agent"],
                "importance": "8-10% JEE weightage. Essential foundation for electrochemistry in Class 12."
            },
            "Hydrogen": {
                "difficulty": "Easy",
                "description": "Hydrogen is the simplest and most abundant element in the universe. This chapter covers hydrogen's unique position in the periodic table (it resembles both alkali metals and halogens), isotopes (protium, deuterium, tritium), preparation and properties of dihydrogen, water structure and properties, heavy water (D₂O), and hydrogen peroxide (H₂O₂). Hydrogen as a future fuel (hydrogen economy) is an important application.",
                "roadmap": ["Position of H in Periodic Table", "Isotopes of Hydrogen", "Preparation & Properties of H₂", "Water & Heavy Water", "Hydrogen Peroxide (H₂O₂)"],
                "topics": ["Electronic configuration 1s¹ — resembles IA and VIIA", "Protium (H), Deuterium (D), Tritium (T)", "Preparation of H₂: electrolysis of water, reaction of metals with acid", "Properties: reducing agent, forms water", "Water: bent shape (V-shaped), bond angle 104.5°, H-bonding", "Hard & Soft water, water treatment", "H₂O₂: mild oxidizing and reducing agent, bleaching", "H₂O₂ structure: non-planar, O-O single bond", "Hydrogen as fuel — fuel cells"],
                "formulas": ["Strength of H₂O₂ in volume strength = 11.2 × Normality", "H₂O₂ volume strength = molarity × 11.2"],
                "applications": ["Hydrogen fuel cells for clean energy", "Heavy water in nuclear reactors", "H₂O₂ in antiseptics and bleaching"],
                "common_mistakes": ["Forgetting H₂O₂ is a covalent compound (not ionic)", "Wrong structure of H₂O₂ (non-planar, not planar)"],
                "importance": "5-7% JEE weightage. Short chapter — factual questions in JEE Main."
            },
            "s-Block Elements (Alkali & Alkaline Earth Metals)": {
                "difficulty": "Medium",
                "description": "s-Block elements are in Groups 1 (alkali metals: Li, Na, K, Rb, Cs, Fr) and Group 2 (alkaline earth metals: Be, Mg, Ca, Sr, Ba, Ra). They have the most reactive metals with characteristic trends in properties. Anomalous behavior of Li (resembles Mg — diagonal relationship) and Be (resembles Al) are important. Compounds like NaOH, Na₂CO₃, NaCl, CaO, CaCO₃, and plaster of paris are industrially vital and frequently tested.",
                "roadmap": ["Group 1 — Alkali Metals", "Group 2 — Alkaline Earth Metals", "Anomalous Behavior of Li & Be", "Important Compounds", "Biological Significance"],
                "topics": ["ns¹ (Group 1) and ns² (Group 2) configurations", "Decreasing IE and EN down the group", "Li anomalous: small size, high charge density, resembles Mg", "Be anomalous: resembles Al (diagonal relationship)", "Flame colors: Li-red, Na-yellow, K-violet, Ca-brick red, Ba-apple green", "NaOH (caustic soda) — chlor-alkali process", "Na₂CO₃ (washing soda), NaHCO₃ (baking soda)", "CaO (quicklime), Ca(OH)₂ (slaked lime), CaCO₃ (limestone)", "Plaster of Paris: 2CaSO₄·H₂O — sets with water"],
                "formulas": ["2Na + 2H₂O → 2NaOH + H₂", "CaO + H₂O → Ca(OH)₂ (slaking)", "2CaSO₄·H₂O + water → CaSO₄·2H₂O (setting of PoP)"],
                "applications": ["Na as heat exchanger in nuclear reactors", "Mg in aircraft alloys", "CaCO₃ in cement production"],
                "common_mistakes": ["Confusing Na₂CO₃ (washing soda) with NaHCO₃ (baking soda)", "Forgetting Be is amphoteric (not basic like other Group 2)"],
                "importance": "8-10% JEE weightage. Important compound properties are frequently tested."
            },
            "p-Block Elements (Groups 13-18)": {
                "difficulty": "Hard",
                "description": "p-Block Elements spans Groups 13 to 18 — the largest section of the periodic table covering diverse elements from reactive non-metals to inert noble gases. This chapter covers boron and aluminium (Group 13), carbon allotropes (Group 14), nitrogen compounds including ammonia and nitric acid (Group 15), sulphuric acid (Group 16), halogens (Group 17), and noble gases (Group 18). Their oxides, oxyacids, and structures are extensively tested in JEE.",
                "roadmap": ["Group 13 — B & Al compounds", "Group 14 — C & Si", "Group 15 — N, P, NH₃, HNO₃", "Group 16 — O, S, H₂SO₄", "Group 17 — Halogens & Group 18 — Noble Gases"],
                "topics": ["Boron: diborane (B₂H₆), boric acid (H₃BO₃), BCl₃", "Aluminium: aluminium chloride, alums (KAl(SO₄)₂·12H₂O)", "Carbon allotropes: diamond, graphite, fullerene", "Si: silicates, zeolites", "Nitrogen: preparation, triple bond, liquid N₂", "NH₃: preparation (Haber process), pyramidal structure, basic", "HNO₃: preparation (Ostwald process), both oxidizer and acid", "Phosphorus allotropes: white, red, black", "PCl₃, PCl₅: structures and hydrolysis", "SO₂, SO₃, H₂SO₄ — contact process", "Ozone (O₃): bent, oxidizing, UV shield", "Halogens: F₂ strongest oxidizer, HF — weak acid", "Interhalogen compounds: ClF₃, BrF₅, IF₇", "Noble gases: Xe compounds — XeF₂, XeF₄, XeF₆"],
                "formulas": ["Haber: N₂ + 3H₂ → 2NH₃ (200 atm, 500°C, Fe catalyst)", "Ostwald: 4NH₃ + 5O₂ → 4NO + 6H₂O (Pt catalyst)", "Contact: 2SO₂ + O₂ → 2SO₃ (V₂O₅ catalyst)"],
                "applications": ["H₂SO₄ — 'king of chemicals' in fertilizer industry", "NH₃ — fertilizers, explosives", "Cl₂ — water purification"],
                "common_mistakes": ["Wrong hybridization for PCl₅ (sp³d) vs PH₃ (sp³)", "Forgetting HF is a weak acid despite F being most electronegative"],
                "importance": "15-18% JEE weightage. Largest chapter — most questions in JEE Main inorganic."
            },
            "Organic Chemistry — Basic Principles & Techniques": {
                "difficulty": "Medium",
                "description": "This chapter is the entry point into organic chemistry — the chemistry of carbon compounds. It covers IUPAC nomenclature, types of isomerism (structural and stereoisomerism), reaction intermediates (carbocations, carbanions, free radicals, carbenes), electronic effects (inductive, mesomeric/resonance, hyperconjugation), and the basic types of organic reactions (substitution, addition, elimination). Mastery here is essential for all subsequent organic chemistry.",
                "roadmap": ["IUPAC Nomenclature", "Structural & Stereoisomerism", "Reaction Intermediates", "Electronic Effects", "Types of Organic Reactions"],
                "topics": ["IUPAC naming of alkanes, alkenes, alkynes, alcohols, aldehydes, acids", "Chain, position, functional group, metamerism isomers", "Geometric (cis-trans) isomerism: condition — different groups on each C", "Optical isomerism: chiral centers, R & S configuration", "Carbocations: stability 3° > 2° > 1°", "Carbanions: stability 1° > 2° > 3°", "Free radicals: stability 3° > 2° > 1°", "Inductive effect: +I and -I groups", "Resonance/mesomeric effect: +M and -M", "Hyperconjugation — stability of alkenes and carbocations", "Electrophilic and nucleophilic substitution, addition, elimination"],
                "formulas": ["+I effect: electron releasing groups", "-I effect: electron withdrawing groups", "+M: lone pairs donate into π system", "-M: withdraw electrons from π system"],
                "applications": ["Drug design (stereoisomers — thalidomide tragedy)", "Reaction mechanism prediction", "Industrial polymer design"],
                "common_mistakes": ["Wrong IUPAC name — not finding longest chain with highest substituents", "Confusing inductive with mesomeric effect", "Wrong stability order of carbocations"],
                "importance": "10-12% JEE weightage. Foundation for entire organic chemistry — master this first."
            },
            "Hydrocarbons": {
                "difficulty": "Hard",
                "description": "Hydrocarbons are organic compounds containing only C and H. They are classified as alkanes (single bonds, CnH₂n₊₂), alkenes (double bond, CnH₂n), alkynes (triple bond, CnH₂n₋₂), and arenes (benzene ring). Each class has characteristic reactions: alkanes undergo free radical substitution, alkenes undergo electrophilic addition (following Markovnikov's rule), alkynes form acetylide salts with sodium, and benzene undergoes electrophilic aromatic substitution.",
                "roadmap": ["Alkanes — Preparation & Free Radical Substitution", "Alkenes — Electrophilic Addition & Markovnikov", "Alkynes — Acidic Character & Addition", "Benzene — Aromaticity & EAS", "Conformations & Cycloalkanes"],
                "topics": ["Alkanes: CnH₂n₊₂, sp³ C, chair & boat cyclohexane", "Preparation: Wurtz reaction, Kolbe's electrolysis, Sabatier-Senderens", "Free radical halogenation: Cl₂/hv, reactivity 3° > 2° > 1°", "Alkenes: sp² C, planar, cis-trans isomers", "Markovnikov's Rule: H adds to more H-bearing C", "Ozonolysis, hydrogenation, hydration, polymerization", "Alkynes: sp C, linear, acetylene", "Lindlar's catalyst: partial hydrogenation → cis alkene", "Na/liquid NH₃: partial hydrogenation → trans alkene", "Acidic H in terminal alkynes, form metal acetylides", "Benzene: 3 double bonds, resonance, Hückel's rule (4n+2 π electrons)", "EAS: nitration, halogenation, sulfonation, Friedel-Crafts alkylation/acylation", "Ortho-para directors (+M groups) vs meta directors (-M groups)"],
                "formulas": ["Hückel: 4n+2 π electrons for aromaticity", "Markovnikov's rule", "Ozonolysis: R₁CH=CHR₂ + O₃ → R₁CHO + R₂CHO"],
                "applications": ["Petrol/diesel (alkanes)", "Ethylene oxide (alkene) in antifreeze", "Benzene in pharmaceuticals and dyes"],
                "common_mistakes": ["Anti-Markovnikov product with HBr + peroxide (free radical)", "Wrong product in ozonolysis - not recognizing ketone vs aldehyde", "Forgetting benzene's EAS (not addition) reactions"],
                "importance": "15-18% JEE weightage. Hydrocarbons reactions are tested heavily every year."
            },
            "Environmental Chemistry": {
                "difficulty": "Easy",
                "description": "Environmental Chemistry deals with chemical processes occurring in the environment and the effects of human activities on air, water, and soil quality. The chapter covers tropospheric and stratospheric pollution, greenhouse effect and global warming, ozone layer depletion mechanism, acid rain, water and soil pollutants, and strategies for environmental protection. Mostly conceptual — important for CBSE boards and JEE Main.",
                "roadmap": ["Atmospheric Pollution", "Greenhouse Effect & Global Warming", "Ozone Layer & Depletion", "Water Pollution", "Soil Pollution & Strategies"],
                "topics": ["Pollutants: primary (CO, SO₂, NO₂, hydrocarbons) and secondary (O₃, PAN)", "CO: toxic, binds haemoglobin, forms HbCO", "Greenhouse gases: CO₂, CH₄, N₂O, CFCs — trap IR radiation", "Global warming: melting ice caps, rising sea levels", "Ozone: O₃ in stratosphere blocks UV-B", "Ozone depletion by CFCs: Cl· radical chain reaction", "Acid rain: H₂SO₄ + HNO₃ from SO₂ + NO₂ in atmosphere", "BOD (Biochemical Oxygen Demand) — water quality measure", "Eutrophication: algal bloom from excess N and P", "Pesticides: DDT, organochlorines — biomagnification"],
                "formulas": ["BOD = dissolved O₂ consumed by bacteria in 5 days at 20°C"],
                "applications": ["CFC phase-out (Montreal Protocol)", "Catalytic converters in cars", "Sewage treatment plants"],
                "common_mistakes": ["Confusing greenhouse effect (thermal blanket) with ozone depletion (UV)", "Not knowing BOD is inversely related to water quality"],
                "importance": "3-5% JEE weightage. Short conceptual chapter — easy marks in JEE Main."
            }
        },

        "Maths": {
            "Trigonometry": {
                "difficulty": "Medium",
                "roadmap": ["Angle Measurement", "Trigonometric Ratios", "Identities & Equations", "Inverse Functions", "Heights & Distances"],
                "topics": ["Radian & Degree Conversion", "Compound Angles", "Multiple & Sub-multiple Angles", "Trigonometric Equations", "Inverse Trig Functions"],
                "formulas": ["sin²θ + cos²θ = 1", "1 + tan²θ = sec²θ", "sin2θ = 2sinθcosθ", "cos2θ = cos²θ - sin²θ"],
                "importance": "Base for calculus. 12-14% board weightage. Essential JEE topic."
            },
            "Sets & Relations": {
                "difficulty": "Easy",
                "roadmap": ["Set Theory Basics", "Operations on Sets", "Relations", "Types of Relations", "Functions"],
                "topics": ["Types of Sets", "Venn Diagrams", "Cartesian Product", "Domain & Range", "Types of Functions"],
                "formulas": ["n(A∪B) = n(A) + n(B) - n(A∩B)", "n(A×B) = n(A) × n(B)"],
                "importance": "Logical foundation. 8% board weightage. Easy scoring."
            },
            "Complex Numbers": {
                "difficulty": "Hard",
                "roadmap": ["Imaginary Numbers", "Algebra of Complex Numbers", "Argand Plane", "Polar Form", "De Moivre's Theorem"],
                "topics": ["i = √(-1)", "Operations on Complex Numbers", "Modulus & Argument", "Euler's Formula", "Roots of Unity"],
                "formulas": ["i² = -1", "|z| = √(a² + b²)", "z = r(cosθ + isinθ)", "(cosθ + isinθ)ⁿ = cosnθ + isinnθ"],
                "importance": "10-12% board weightage. High-value JEE problems."
            },
            "Permutations & Combinations": {
                "difficulty": "Medium",
                "roadmap": ["Fundamental Principle", "Permutations", "Combinations", "Circular Arrangements", "Restricted Cases"],
                "topics": ["Multiplication & Addition Principle", "nPr", "nCr", "Circular Permutations", "Selection with Constraints"],
                "formulas": ["nPr = n!/(n-r)!", "nCr = n!/[r!(n-r)!]", "nCr = nCn-r"],
                "importance": "10% board weightage. Probability base. Appears in JEE."
            },
            "Binomial Theorem": {
                "difficulty": "Medium",
                "roadmap": ["Expansion Formula", "General & Middle Terms", "Binomial Coefficients", "Applications"],
                "topics": ["(a+b)ⁿ Expansion", "Tr+1 Term", "Greatest Coefficient", "Pascal's Triangle"],
                "formulas": ["(a+b)ⁿ = Σ nCr·aⁿ⁻ʳ·bʳ", "Tr+1 = nCr·aⁿ⁻ʳ·bʳ"],
                "importance": "8-10% board weightage. Simple mark-scoring chapter."
            },
            "Sequences & Series": {
                "difficulty": "Hard",
                "roadmap": ["AP - Arithmetic Progression", "GP - Geometric Progression", "HP - Harmonic Progression", "Special Series", "AM, GM, HM"],
                "topics": ["nth Term Formulas", "Sum of n Terms", "Infinite GP", "Arithmetic Mean", "Geometric Mean Inequality"],
                "formulas": ["Tn = a + (n-1)d", "Sn = n/2[2a + (n-1)d]", "Tn = arⁿ⁻¹", "S∞ = a/(1-r) for |r|<1", "GM² = AM × HM"],
                "importance": "12% board weightage. Crucial for series convergence in calculus."
            },
            "Straight Lines": {
                "difficulty": "Medium",
                "roadmap": ["Distance & Section Formula", "Slope of Line", "Various Forms of Line", "Angle Between Lines", "Distance from Point to Line"],
                "topics": ["Slope-Intercept Form", "Point-Slope Form", "Two-Point Form", "Parallel & Perpendicular Lines", "Foot of Perpendicular"],
                "formulas": ["m = (y₂-y₁)/(x₂-x₁)", "y = mx + c", "ax + by + c = 0", "d = |ax₁+by₁+c|/√(a²+b²)"],
                "importance": "10% board weightage. Geometry foundation for JEE."
            },
            "Conic Sections": {
                "difficulty": "Hard",
                "roadmap": ["Circle", "Parabola", "Ellipse", "Hyperbola", "General Equation"],
                "topics": ["Standard Forms", "Focus & Directrix", "Eccentricity", "Latus Rectum", "Tangent & Normal"],
                "formulas": ["x² + y² = r²", "y² = 4ax", "x²/a² + y²/b² = 1", "x²/a² - y²/b² = 1"],
                "importance": "15% board weightage. Maximum JEE questions from coordinate geometry."
            },
            "Limits & Derivatives": {
                "difficulty": "Hard",
                "roadmap": ["Concept of Limits", "Evaluation Techniques", "Continuity", "Derivative Definition", "Differentiation Rules"],
                "topics": ["Left & Right Limits", "L'Hospital's Rule", "Continuous Functions", "First Principles", "Product, Quotient, Chain Rule"],
                "formulas": ["lim(x→0) sinx/x = 1", "d/dx(xⁿ) = nxⁿ⁻¹", "d/dx(uv) = u'v + uv'"],
                "importance": "Foundation for Class 12. 12% board weightage. Core calculus concept."
            }
        }
    },
    "12": {
        "Physics": {
            "Electric Charges & Fields": {
                "difficulty": "Hard",
                "description": "Electric Charges & Fields is the first chapter of Class 12 Physics and the foundation of electrostatics. It introduces the concept of electric charge, Coulomb's law for force between charges, and the electric field — a vector field that represents the force per unit charge at any point in space. Gauss's Law is a powerful theorem that relates the electric flux through a closed surface to the enclosed charge, simplifying field calculations for symmetric configurations like spheres, cylinders, and planes.",
                "roadmap": ["Properties of Electric Charge", "Coulomb's Law", "Electric Field & Field Lines", "Electric Dipole", "Gauss's Law & Applications"],
                "topics": ["Charge Quantization: q = ne", "Superposition Principle", "Electric Field E = F/q", "Field due to dipole (axial & equatorial)", "Electric Flux Φ = E·A = E A cosθ", "Gauss's Law: ΦE = q_enc/ε₀", "Field due to infinite line charge, sheet, sphere"],
                "formulas": ["F = kq₁q₂/r² = q₁q₂/4πε₀r²", "E = kq/r²", "E_dipole(axial) = 2kp/r³", "E_dipole(equatorial) = kp/r³", "Φ = q/ε₀ (Gauss's Law)"],
                "laws": ["Coulomb's Law", "Gauss's Law"],
                "applications": ["Lightning rods", "Xerography (photocopying)", "Electrostatic precipitators"],
                "common_mistakes": ["Forgetting Gaussian surface must be closed for Gauss's Law", "Direction errors in dipole field"],
                "importance": "15-18% JEE weightage. Gauss's Law problems and dipole questions are JEE staples."
            },
            "Electrostatic Potential & Capacitance": {
                "difficulty": "Hard",
                "description": "Electrostatic potential is the work done per unit charge in bringing a test charge from infinity to a given point in an electric field. It is a scalar quantity, making calculations simpler than vector field problems. Capacitors store energy in electric fields and are critical components in circuits. This chapter covers potential due to charges and dipoles, equipotential surfaces, dielectrics, and combinations of capacitors — all heavily tested in JEE.",
                "roadmap": ["Electric Potential & Relation to E", "Potential due to a Dipole", "Equipotential Surfaces", "Capacitors & Capacitance", "Dielectrics & Energy Stored"],
                "topics": ["V = kq/r (point charge)", "V due to dipole: V = kp cosθ/r²", "Relation: E = −dV/dr", "Capacitance C = Q/V", "Parallel plate: C = ε₀A/d", "Series: 1/C = 1/C₁ + 1/C₂", "Parallel: C = C₁ + C₂", "Energy U = ½CV² = Q²/2C", "Dielectric constant K, C' = KC"],
                "formulas": ["V = kq/r", "E = −∇V", "C = ε₀A/d", "U = ½CV²", "U = ε₀E²/2 (energy density)"],
                "applications": ["Camera flash (capacitor discharge)", "Power factor correction", "Touchscreens"],
                "common_mistakes": ["Confusing potential (scalar) with field (vector)", "Mixing up series/parallel capacitor formulas"],
                "importance": "High JEE weightage. Capacitor combinations and energy calculations appear every year."
            },
            "Current Electricity": {
                "difficulty": "Medium",
                "description": "Current Electricity deals with the physics of electric current — the flow of charge through conductors. Ohm's Law, resistance, resistivity, and their temperature dependence form the core. Kirchhoff's Laws (KVL and KCL) are the tools for solving complex circuit networks. The Wheatstone bridge principle is used in the Meter Bridge and Potentiometer — instruments that appear in every JEE paper. Understanding EMF, internal resistance, and maximum power transfer is essential for circuit analysis.",
                "roadmap": ["Electric Current & Drift Velocity", "Ohm's Law & Resistance", "Temperature Dependence", "Kirchhoff's Laws", "Wheatstone Bridge & Potentiometer"],
                "topics": ["Drift velocity vd = I/nAe", "Resistivity ρ, R = ρl/A", "Colour code of resistors", "KVL: ΣV = 0 in loop", "KCL: ΣI = 0 at node", "Wheatstone condition: P/Q = R/S", "EMF & internal resistance: V = E − Ir", "Maximum power transfer: r = R_external", "Potentiometer for comparing EMFs"],
                "formulas": ["I = nAevd", "V = IR", "R = ρl/A", "Rₜ = R₀(1 + αT)", "P = VI = I²R = V²/R", "η_max = 50% (max power transfer)"],
                "laws": ["Ohm's Law", "Kirchhoff's Current Law (KCL)", "Kirchhoff's Voltage Law (KVL)"],
                "applications": ["Home wiring circuits", "Battery charging systems", "Resistance thermometers"],
                "common_mistakes": ["Wrong sign in KVL loop equations", "Forgetting internal resistance reduces terminal voltage"],
                "importance": "15% JEE weightage. Kirchhoff's Laws and Wheatstone bridge are annual JEE questions."
            },
            "Moving Charges & Magnetism": {
                "difficulty": "Hard",
                "description": "This chapter explores how moving electric charges create magnetic fields and how magnetic fields exert forces on moving charges. The Biot-Savart Law gives the magnetic field due to a current element, while Ampere's Law provides a simpler method for symmetric configurations. The Lorentz force on a charged particle in a magnetic field causes circular or helical motion — the basis of cyclotrons and mass spectrometers. Galvanometers, motors, and current-carrying conductors in fields are key JEE topics.",
                "roadmap": ["Biot-Savart Law", "Ampere's Circuital Law", "Force on a Charge & Current in B", "Torque on Current Loop", "Moving Coil Galvanometer"],
                "topics": ["Biot-Savart: dB = μ₀Idl sinθ/4πr²", "Field due to straight wire: B = μ₀I/2πr", "Field at center of circular loop: B = μ₀I/2R", "Solenoid: B = μ₀nI", "Lorentz force: F = q(v×B)", "Cyclotron frequency: f = qB/2πm", "Torque: τ = NIAB sinθ", "Galvanometer: current sensitivity = NBA/k"],
                "formulas": ["B = μ₀I/2πr (wire)", "B = μ₀nI (solenoid)", "F = qvB sinθ", "F = BIl sinθ (on conductor)", "r = mv/qB (circular motion)"],
                "laws": ["Biot-Savart Law", "Ampere's Circuital Law"],
                "applications": ["Electric motors", "Cyclotron particle accelerator", "MRI machines"],
                "common_mistakes": ["Confusing direction of magnetic force (use Fleming's left-hand rule)", "Forgetting sin θ factor"],
                "importance": "15% JEE weightage. Numericals on field calculations and force on conductors are common."
            },
            "Magnetism & Matter": {
                "difficulty": "Medium",
                "description": "Magnetism & Matter studies the magnetic properties of materials — why some are attracted to magnets (ferromagnetics), some are weakly attracted (paramagnetics), and some are weakly repelled (diamagnetics). The Earth's magnetic field and its elements (declination, dip, horizontal component) are important concepts. Magnetic dipoles, their moment, potential energy in a field, and the hysteresis loop of ferromagnetic materials are also covered in this chapter.",
                "roadmap": ["Bar Magnet as Magnetic Dipole", "Earth's Magnetism", "Magnetic Properties of Materials", "Hysteresis Loop", "Permanent Magnets & Electromagnets"],
                "topics": ["Magnetic moment m = IA", "Field of bar magnet (axial & equatorial)", "Earth's magnetic elements: Declination, Dip, H", "Dia, Para, Ferromagnetism", "Curie Temperature", "B-H curve & Hysteresis loss", "Retentivity & Coercivity"],
                "formulas": ["B_axial = μ₀/4π × 2M/r³", "B_equatorial = μ₀/4π × M/r³", "U = −m·B = −mB cosθ", "χ_m = M/H (susceptibility)"],
                "applications": ["Permanent magnets in speakers", "Transformers (low hysteresis loss needed)", "Magnetic recording media"],
                "common_mistakes": ["Confusing B (magnetic flux density) with H (magnetic intensity)", "Mixing up dia/para/ferromagnetic properties"],
                "importance": "5-8% JEE weightage. Earth's magnetism and classification problems appear in JEE Main."
            },
            "Electromagnetic Induction": {
                "difficulty": "Hard",
                "description": "Electromagnetic Induction (EMI) is one of the most important discoveries in physics — the basis of electricity generation worldwide. Faraday's Law states that changing magnetic flux induces an EMF. Lenz's Law gives its direction (opposing the change). Motional EMF in a conductor moving in a magnetic field, self-inductance (L), mutual inductance (M), and the energy stored in an inductor are all key concepts. AC generators, transformers, and eddy currents are critical applications tested every year in JEE.",
                "roadmap": ["Magnetic Flux & Faraday's Law", "Lenz's Law", "Motional EMF", "Self & Mutual Inductance", "AC Generator & Transformer"],
                "topics": ["Φ = BA cosθ", "ε = −dΦ/dt (Faraday)", "Lenz's Law — direction of induced current", "Motional EMF ε = Blv", "Self inductance L = Φ/I, ε = −L(dI/dt)", "Mutual inductance M = Φ₂/I₁", "Energy in inductor: U = ½LI²", "AC generator: ε = NBAω sinωt", "Transformer: Vs/Vp = Ns/Np"],
                "formulas": ["ε = −NdΦ/dt", "ε = Blv", "L_solenoid = μ₀n²Al", "U = ½LI²", "Vs/Vp = Ns/Np = Ip/Is"],
                "laws": ["Faraday's Law of Electromagnetic Induction", "Lenz's Law"],
                "applications": ["Electric generators & alternators", "Transformers in power grid", "Induction cooktops (eddy currents)"],
                "common_mistakes": ["Sign error in Lenz's law direction", "Confusing L (self) with M (mutual) inductance"],
                "importance": "18-20% JEE weightage. Highest scoring chapter in Class 12 Physics for JEE."
            },
            "Alternating Current": {
                "difficulty": "Hard",
                "description": "Alternating Current (AC) circuits are the backbone of modern electrical infrastructure. Unlike DC, AC current and voltage vary sinusoidally with time. This chapter introduces RMS values, impedance of R, L, C elements, phasors, resonance in LC and RLC circuits, and the power factor. The LCR series circuit at resonance is a particularly important concept because it has the highest impedance selectivity and is tested every year. Transformers and their efficiency are also covered here.",
                "roadmap": ["AC Voltage & Current — RMS Values", "Phasors & Phase Relationships", "Reactance of L and C", "LCR Series Circuit & Resonance", "Power in AC Circuits"],
                "topics": ["V_rms = V₀/√2, I_rms = I₀/√2", "Reactance: X_L = ωL, X_C = 1/ωC", "Impedance Z = √(R²+(XL−XC)²)", "Phase angle: tan φ = (XL−XC)/R", "Resonance: ω₀ = 1/√LC, Z = R (minimum)", "Quality factor Q = ω₀L/R", "Power P = V_rms I_rms cos φ", "Power factor cos φ = R/Z"],
                "formulas": ["Z = √(R² + (XL−XC)²)", "ω₀ = 1/√(LC)", "Q = ω₀L/R = 1/ω₀CR", "P_avg = V_rms I_rms cosφ", "V_rms = V₀/√2"],
                "applications": ["Radio tuning (resonance)", "Power factor correction in industries", "AC power transmission"],
                "common_mistakes": ["Confusing peak and RMS values", "Forgetting resonance condition XL = XC"],
                "importance": "12-15% JEE weightage. Resonance, impedance, and power factor appear every year in JEE Main."
            },
            "Electromagnetic Waves": {
                "difficulty": "Easy",
                "description": "Electromagnetic Waves are transverse waves consisting of oscillating electric and magnetic fields perpendicular to each other and to the direction of propagation. They do not require a medium and travel at the speed of light c = 3×10⁸ m/s in vacuum. Maxwell's equations predicted their existence. The electromagnetic spectrum ranges from radio waves (lowest frequency) to gamma rays (highest frequency). The seven types of EM waves and their sources and uses are commonly tested in JEE.",
                "roadmap": ["Maxwell's Equations & Displacement Current", "Properties of EM Waves", "Electromagnetic Spectrum", "Speed of EM Waves"],
                "topics": ["Displacement current I_d = ε₀ dΦE/dt", "E & B perpendicular to each other & to propagation", "c = E₀/B₀ = 1/√(μ₀ε₀)", "γ-rays, X-rays, UV, Visible, IR, Microwaves, Radio waves", "Sources & uses of each type", "Wave intensity & energy density"],
                "formulas": ["c = E/B = 1/√(μ₀ε₀) ≈ 3×10⁸ m/s", "u = ½ε₀E² + B²/2μ₀ (energy density)", "I_d = ε₀ dΦE/dt"],
                "laws": ["Maxwell's Equations (conceptual understanding)"],
                "applications": ["Radio/TV broadcasting", "X-ray medical imaging", "Microwave cooking"],
                "common_mistakes": ["Confusing frequency ordering in EM spectrum", "Forgetting E and B are in phase in EM waves"],
                "importance": "5-8% JEE weightage. Easy chapter — mostly conceptual. Know the EM spectrum order."
            },
            "Ray Optics & Optical Instruments": {
                "difficulty": "Medium",
                "description": "Ray Optics treats light as travelling in straight lines (rays) and studies reflection and refraction. Mirrors and lenses form images using well-defined formulas. Total Internal Reflection — the basis of optical fibres — occurs when light goes from a denser to a rarer medium beyond the critical angle. Optical instruments — microscopes, telescopes, and the human eye — are built from basic lens/mirror combinations. Prisms and dispersion of light are also important JEE topics.",
                "roadmap": ["Reflection — Mirrors & Image Formation", "Refraction — Snell's Law & TIR", "Prism & Dispersion", "Lenses — Convex, Concave", "Optical Instruments"],
                "topics": ["Mirror formula: 1/f = 1/v + 1/u", "Magnification m = −v/u", "Snell's Law: μ₁ sinθ₁ = μ₂ sinθ₂", "TIR: sin C = 1/μ (critical angle)", "Lens formula: 1/f = 1/v − 1/u", "Lens maker's equation: 1/f = (μ−1)(1/R₁ − 1/R₂)", "Power P = 1/f (in metres)", "Prism: deviation δ = A(μ−1)", "Microscope & telescope magnification"],
                "formulas": ["1/v − 1/u = 1/f (lens)", "1/v + 1/u = 1/f (mirror)", "μ = c/v = sin i/sin r", "P = 1/f, P = P₁ + P₂ (combined)"],
                "applications": ["Optical fibres (TIR)", "Cameras & projectors", "Corrective lenses (spectacles)"],
                "common_mistakes": ["Mixing sign conventions for mirrors vs lenses", "Using mirror formula for lens or vice versa"],
                "importance": "12-15% JEE weightage. Lens/mirror combinations and TIR problems are very common."
            },
            "Wave Optics": {
                "difficulty": "Hard",
                "description": "Wave Optics treats light as a wave and explains phenomena that Ray Optics cannot — interference, diffraction, and polarisation. Young's Double Slit Experiment (YDSE) demonstrates interference and is one of the most important experiments in physics. Single-slit diffraction shows central maxima and minima. Polarisation by reflection (Brewster's angle) and by polaroids has applications in LCD screens and sunglasses. Coherence, wavefronts, and Huygens' principle underlie all wave optics.",
                "roadmap": ["Huygens' Principle & Wavefronts", "Young's Double Slit — Interference", "Single Slit Diffraction", "Resolving Power", "Polarisation"],
                "topics": ["Coherent sources for sustained interference", "Path difference Δ = d sinθ ≈ yd/D", "Fringe width β = λD/d", "Constructive: Δ = nλ, Destructive: Δ = (n+½)λ", "Single slit: minima at a sinθ = nλ", "Diffraction limit", "Brewster's angle: tan θ_B = μ", "Malus's Law: I = I₀ cos²θ"],
                "formulas": ["β = λD/d", "y_n = nλD/d", "tan θ_B = μ (Brewster)", "I = I₀ cos²θ (Malus)"],
                "applications": ["Anti-reflection coating", "Holography", "Polaroid sunglasses"],
                "common_mistakes": ["Confusing path difference with phase difference", "Forgetting fringe width is independent of order n"],
                "importance": "12-15% JEE weightage. YDSE is the most tested optics topic. Know fringe width formula by heart."
            },
            "Dual Nature of Radiation & Matter": {
                "difficulty": "Hard",
                "description": "Dual Nature of Radiation & Matter established that light behaves as both a wave and a particle (photon), and conversely, matter (electrons) can behave as waves. Einstein's explanation of the Photoelectric Effect — for which he won the Nobel Prize — showed that light comes in quanta of energy hν. de Broglie hypothesised that moving particles have a wavelength λ = h/p, confirmed by electron diffraction. This chapter is the bridge between classical and quantum physics.",
                "roadmap": ["Photoelectric Effect — Einstein's Theory", "Photon & its Properties", "Threshold Frequency & Work Function", "de Broglie Hypothesis", "Davisson-Germer Experiment"],
                "topics": ["Photon energy E = hν = hc/λ", "Photoelectric equation: KE_max = hν − φ", "Threshold frequency ν₀ = φ/h", "Stopping potential eV₀ = KE_max", "de Broglie wavelength λ = h/p = h/mv", "λ for accelerated electron: λ = h/√(2meV)", "Wave-particle duality"],
                "formulas": ["E = hν", "KE_max = hν − φ", "λ = h/mv (de Broglie)", "λ = h/√(2meV) (electron in field V)"],
                "laws": ["Einstein's Photoelectric Law"],
                "applications": ["Solar cells", "Photomultiplier tubes", "Electron microscopes (de Broglie)"],
                "common_mistakes": ["Confusing frequency with wavelength for threshold", "Forgetting stopping potential ≠ KE directly"],
                "importance": "10-12% JEE weightage. Photoelectric effect and de Broglie wavelength appear every year."
            },
            "Atoms": {
                "difficulty": "Medium",
                "description": "The Atoms chapter traces the development of atomic models culminating in Bohr's model of the hydrogen atom. Rutherford's gold foil experiment proved the nuclear model — an atom with a tiny, dense, positive nucleus. Bohr quantised the electron's angular momentum and derived energy levels, explaining hydrogen's line spectra. The transitions between energy levels correspond to emission or absorption of photons of specific wavelengths — giving rise to spectral series like Lyman, Balmer, Paschen.",
                "roadmap": ["Thomson's & Rutherford's Models", "Bohr's Postulates & Energy Levels", "Hydrogen Spectrum — Spectral Series", "Limitations of Bohr's Model"],
                "topics": ["Rutherford's nuclear model & scattering experiment", "Bohr's postulates: mvr = nh/2π", "Energy of nth orbit: Eₙ = −13.6/n² eV", "Radius: rₙ = n²a₀ (a₀ = 0.529 Å)", "Spectral series: Lyman(UV), Balmer(visible), Paschen(IR)", "1/λ = R(1/n₁² − 1/n₂²)"],
                "formulas": ["Eₙ = −13.6/n² eV", "rₙ = 0.529 n² Å", "1/λ = R_H(1/n₁² − 1/n₂²)", "v = 2.18×10⁶/n m/s"],
                "applications": ["Atomic clocks", "Laser (stimulated emission)", "Hydrogen fuel cells"],
                "common_mistakes": ["Confusing emission (photon released) with absorption (photon absorbed)", "Wrong series — Lyman ends at n=1, Balmer at n=2"],
                "importance": "8-10% JEE weightage. Spectral series and Bohr energy calculations are frequent questions."
            },
            "Nuclei": {
                "difficulty": "Hard",
                "description": "The Nuclei chapter dives into the atomic nucleus — a tiny but incredibly dense core containing protons and neutrons (nucleons). Nuclear forces holding the nucleus together are the strongest forces in nature. Mass defect and binding energy explain why nuclei are stable. Radioactive nuclei decay by alpha (α), beta (β), or gamma (γ) emission following a precise exponential decay law. Nuclear fission (splitting heavy nuclei) and fusion (joining light nuclei) release enormous energy — the basis of nuclear reactors and the Sun.",
                "roadmap": ["Nucleus — Size, Mass, Composition", "Mass Defect & Binding Energy", "Radioactivity — α, β, γ Decay", "Radioactive Decay Law & Half-Life", "Nuclear Fission & Fusion"],
                "topics": ["Nuclear density ≈ 2.3×10¹⁷ kg/m³", "Mass defect Δm, BE = Δmc²", "BE per nucleon curve", "Activity A = λN = A₀e^(−λt)", "Half life T½ = 0.693/λ", "N = N₀(½)^(t/T½)", "Conservation laws in decay", "Q-value of nuclear reaction", "Chain reaction & critical mass"],
                "formulas": ["BE = (Zm_p + Nm_n − M)c²", "N = N₀ e^(−λt)", "T½ = 0.693/λ", "A = λN", "Q = (m_reactants − m_products)c²"],
                "laws": ["Radioactive Decay Law (exponential)"],
                "applications": ["Nuclear power plants (fission)", "Carbon dating (radioactive decay)", "Hydrogen bomb (fusion)"],
                "common_mistakes": ["Confusing half-life with mean life (τ = T½/0.693)", "Forgetting to use 931.5 MeV/amu for mass-energy"],
                "importance": "10-12% JEE weightage. Binding energy curve and decay calculations appear every year."
            },
            "Semiconductor Electronics": {
                "difficulty": "Medium",
                "description": "Semiconductor Electronics is the chapter that explains the devices powering our digital world. Semiconductors have conductivity between conductors and insulators. Doping creates n-type (extra electrons) and p-type (extra holes) semiconductors. The p-n junction diode allows current in one direction (forward bias) and blocks in reverse — the basis of rectifiers. Transistors act as amplifiers and switches. Logic gates (AND, OR, NOT, NAND, NOR, XOR) form the foundation of all digital electronics.",
                "roadmap": ["Energy Bands — Conductor, Insulator, Semiconductor", "p-type & n-type Semiconductors", "p-n Junction Diode & Rectification", "Transistor — CE Configuration", "Logic Gates & Boolean Algebra"],
                "topics": ["Intrinsic & extrinsic semiconductors", "Fermi level concept", "Depletion layer & barrier potential", "Forward & reverse biased diode I-V characteristics", "Half-wave and full-wave rectifier", "Zener diode as voltage regulator", "BJT in CE mode: current gain β = Ic/Ib", "NAND & NOR as universal gates", "Truth tables of all basic gates"],
                "formulas": ["I = I₀(e^(eV/kT) − 1) (diode)", "β = Ic/Ib", "α = Ic/Ie", "β = α/(1−α)", "A_V = β(Rc/Ri) (CE amplifier)"],
                "applications": ["Mobile phones & computers (transistors)", "Solar cells (p-n junction)", "LED lighting"],
                "common_mistakes": ["Confusing α and β for transistor", "Forgetting NAND = NOT AND (outputs are inverted)"],
                "importance": "10% JEE weightage. Logic gates and diode rectification are easy marks in JEE Main."
            },
            "Communication Systems": {
                "difficulty": "Easy",
                "description": "Communication Systems is the most conceptual and least mathematical chapter in Class 12 Physics. It covers how information is transmitted from one place to another using electromagnetic signals. The key processes are modulation (combining information with a carrier wave) and demodulation. Amplitude modulation (AM), frequency modulation (FM), and digital communication are the main types. The chapter also covers bandwidth, signal attenuation, modes of propagation, and the internet/mobile network ecosystem.",
                "roadmap": ["Elements of Communication System", "Signal Bandwidth & Channel", "Modulation — AM & FM", "Signal Propagation Modes", "Internet, Mobile & Satellite Communication"],
                "topics": ["Source, Transmitter, Channel, Receiver, Sink", "Bandwidth of AM signal = 2f_m", "Modulation index m = Am/Ac", "Ground wave, Sky wave, Space wave propagation", "Critical frequency fc = 9√Nmax (ionosphere)", "Maximum usable frequency (MUF)", "Satellite communication bands"],
                "formulas": ["Bandwidth = 2f_m (AM)", "m = Am/Ac (modulation index)", "Range of TV tower d = √(2Rh)"],
                "applications": ["AM/FM radio broadcasting", "Satellite TV and GPS", "Mobile communication (4G/5G)"],
                "common_mistakes": ["Confusing AM bandwidth formula", "Mixing up ground wave vs sky wave conditions"],
                "importance": "5-8% JEE weightage. Theory-based easy marks. Know AM bandwidth and modulation index."
            }
        },
        "Chemistry": {
            "Electrochemistry": {
                "difficulty": "Hard",
                "description": "Electrochemistry connects chemistry and electricity. Galvanic (voltaic) cells convert chemical energy to electrical energy spontaneously; electrolytic cells use electrical energy to drive non-spontaneous chemical reactions (electrolysis). The Nernst Equation adjusts standard cell potential for non-standard conditions. Conductance, Kohlrausch's Law, and batteries (dry cell, lead-acid, Li-ion, fuel cells) are all key topics tested in JEE annually.",
                "roadmap": ["Redox Reactions Review", "Galvanic & Electrolytic Cells", "Standard Electrode Potential & EMF", "Nernst Equation", "Conductance & Batteries"],
                "topics": ["Galvanic cell: anode (oxidation, -) and cathode (reduction, +)", "Salt bridge: maintains electrical neutrality", "Standard hydrogen electrode (SHE): E° = 0 V", "E°cell = E°cathode - E°anode", "Electrochemical series — relative reactivity", "Nernst equation: E = E° - (RT/nF)ln Q = E° - (0.0591/n)log Q at 25°C", "ΔG = -nFE; ΔG° = -nFE° = -RT ln K", "Faraday's 1st law: m = ZIt = (M/nF)It", "Faraday's 2nd law: m₁/m₂ = E₁/E₂ (equivalent weights)", "Conductance G = 1/R", "Specific conductance κ = G × l/A", "Molar conductance Λm = κ/c × 1000", "Kohlrausch's Law: Λ°m = Σλ°(cations) + Σλ°(anions)", "Dry cell: Zn-MnO₂", "Lead storage battery: Pb-PbO₂-H₂SO₄", "Fuel cells: H₂-O₂ cell (continuous supply)", "Corrosion: Fe oxidized to Fe²⁺ → rust (Fe₂O₃·xH₂O)"],
                "formulas": ["E°cell = E°cathode - E°anode", "E = E° - (0.0591/n)log Q (Nernst)", "ΔG = -nFE", "m = (M × I × t)/(n × F) (Faraday)", "Λm = κ × 1000/c"],
                "applications": ["Lead-acid batteries in vehicles", "Electroplating (silver, gold)", "Fuel cells in spacecraft"],
                "common_mistakes": ["Confusing anode/cathode sign in different cell types", "Forgetting temperature in Nernst equation", "Wrong formula for Faraday's Law (using n for moles, not equivalents)"],
                "importance": "15-18% JEE weightage. Nernst equation and electrolysis problems appear every year."
            },
            "Chemical Kinetics": {
                "difficulty": "Hard",
                "description": "Chemical Kinetics studies the rate (speed) of chemical reactions and the factors that affect it. The rate law relates reaction rate to reactant concentrations. Integrated rate laws allow calculation of concentration at any time. Arrhenius Equation connects rate constant to temperature and activation energy. Collision theory and transition state theory explain the molecular basis of reactions. Catalysts speed up reactions by providing an alternative pathway with lower activation energy.",
                "roadmap": ["Rate of Reaction & Factors Affecting It", "Rate Law & Order of Reaction", "Integrated Rate Equations", "Temperature Dependence — Arrhenius", "Collision Theory & Catalysis"],
                "topics": ["Rate = -d[A]/dt = +d[B]/dt (stoichiometry)", "Rate law: Rate = k[A]^m[B]^n", "Order = m + n (from experiment, not stoichiometry)", "Molecularity: unimolecular, bimolecular, termolecular", "Zero order: [A] = [A]₀ - kt; t½ = [A]₀/2k", "First order: ln[A] = ln[A]₀ - kt; t½ = 0.693/k", "Units of k: zero order (mol L⁻¹ s⁻¹), first order (s⁻¹)", "Arrhenius: k = Ae^(-Ea/RT)", "log(k₂/k₁) = Ea/2.303R × (T₂-T₁)/(T₁T₂)", "Catalysis: homogeneous, heterogeneous, enzyme (Michaelis-Menten)", "Activation energy Ea", "Pseudo-first order reactions"],
                "formulas": ["Rate = k[A]^m[B]^n", "t½ = 0.693/k (first order)", "k = Ae^(-Ea/RT) (Arrhenius)", "log(k₂/k₁) = Ea/2.303R × (T₂-T₁)/(T₁T₂)"],
                "applications": ["Refrigeration (slowing food spoilage)", "Industrial catalyst design", "Drug metabolism rates in pharmacology"],
                "common_mistakes": ["Determining order from balanced equation (order is from experiment!)", "Wrong units of k for different orders", "Confusing t½ formulas for zero, first, second order"],
                "importance": "12-15% JEE weightage. Graph-based questions (ln[A] vs t, [A] vs t) are very common."
            },
            "Surface Chemistry": {
                "difficulty": "Medium",
                "description": "Surface Chemistry deals with phenomena occurring at surfaces and interfaces. Adsorption explains how substances accumulate at surfaces (physosorption vs chemisorption). Colloids are heterogeneous mixtures with particle sizes 1-1000 nm that show unique properties like the Tyndall effect and Brownian motion. Emulsions, catalysis, and enzyme action are practical applications. This chapter is mostly conceptual and theory-based.",
                "roadmap": ["Adsorption — Physosorption vs Chemisorption", "Freundlich & Langmuir Adsorption Isotherms", "Catalysis", "Colloidal State — Types & Properties", "Emulsions & Applications"],
                "topics": ["Adsorption: surface phenomenon, Δ H < 0", "Physosorption: weak van der Waals, reversible, multilayer", "Chemisorption: strong chemical bond, irreversible, monolayer", "Freundlich adsorption isotherm: x/m = kp^(1/n)", "Catalysis: homogeneous (same phase), heterogeneous (different phase)", "Colloidal particle size: 1-1000 nm", "Tyndall effect: scattering of light by colloidal particles", "Brownian motion: random movement due to molecular bombardment", "Electrophoresis: movement of colloidal particles in electric field", "Coagulation/flocculation: destabilization of colloid", "Hardy-Schulze rule: higher valency of coagulating ion, faster coagulation", "Types of colloids: lyophilic (stable) and lyophobic (unstable)", "Micelles: soap/detergent in water above CMC", "Emulsions: oil-in-water and water-in-oil"],
                "formulas": ["x/m = kp^(1/n) (Freundlich)", "x/m = ap/(1+bp) (Langmuir)"],
                "applications": ["Catalytic converters (heterogeneous catalysis)", "Dialysis for kidney patients", "Food emulsions (mayonnaise, milk)"],
                "common_mistakes": ["Confusing lyophilic and lyophobic colloids", "Forgetting Tyndall effect is shown by colloids NOT true solutions", "Wrong Hardy-Schulze rule application"],
                "importance": "8-10% JEE weightage. Mostly conceptual — 2-3 one-liner questions in JEE Main."
            },
            "d & f Block Elements": {
                "difficulty": "Hard",
                "description": "d-Block elements (transition metals, Groups 3-12) have partially filled d-orbitals giving them unique properties: variable oxidation states, colored compounds, magnetic properties, and catalytic activity. f-Block elements (lanthanoids and actinoids) have 4f and 5f orbitals being filled. K₂Cr₂O₇ and KMnO₄ are the most important oxidizing agents tested in JEE. Lanthanoid contraction and its consequences are crucial concepts.",
                "roadmap": ["Electronic Configuration — (n-1)d¹⁻¹⁰ns¹⁻²", "Properties of Transition Metals", "Important Compounds: K₂Cr₂O₇ & KMnO₄", "Lanthanoids & Lanthanoid Contraction", "Actinoids & Comparison"],
                "topics": ["d-block: (n-1)d¹⁻¹⁰ns¹⁻² variable oxidation states", "Variable valency: multiple d-electrons can be involved in bonding", "Color due to d-d transitions in unpaired electrons", "Magnetic moment: μ = √n(n+2) BM (n = unpaired e⁻)", "KMnO₄: Mn goes from +7 to +2 (acid), +4 (neutral), +6 (base)", "K₂Cr₂O₇: Cr goes from +6 to +3 in acidic reactions", "Cr₂O₇²⁻ + 14H⁺ + 6e⁻ → 2Cr³⁺ + 7H₂O", "Lanthanoids: La to Lu, 4f filling, Ce-Eu and Gd-Yb", "Lanthanoid contraction: steady decrease in atomic radius due to poor 4f shielding", "Consequence: Zr ≈ Hf (same radius), similar chemistry", "Actinoids: Ac to Lr, 5f filling, all radioactive"],
                "formulas": ["μ = √n(n+2) BM (magnetic moment)", "E = hν = hc/λ (energy of light absorbed, gives color)"],
                "applications": ["Fe catalyst in Haber process", "Pt/Pd in catalytic converters", "Ti alloys in aerospace (titanium)"],
                "common_mistakes": ["Wrong electronic config of Cr (3d⁵4s¹) and Cu (3d¹⁰4s¹)", "Forgetting KMnO₄ product changes with medium (acidic/neutral/basic)", "Wrong magnetic moment calculation"],
                "importance": "12-15% JEE weightage. KMnO₄/K₂Cr₂O₇ reactions and magnetic moment appear every year."
            },
            "Coordination Compounds": {
                "difficulty": "Hard",
                "description": "Coordination Compounds are compounds where a central metal atom/ion is bonded to surrounding ligands via coordinate (dative) bonds. Werner's Theory first explained their structure. IUPAC nomenclature, types of isomerism (geometrical, optical, ionization, linkage), Valence Bond Theory (VBT) for hybridization, and Crystal Field Theory (CFT) for color and magnetism are all tested extensively in JEE.",
                "roadmap": ["Werner's Theory & Terminology", "IUPAC Nomenclature", "Isomerism in Coordination Compounds", "Valence Bond Theory (VBT)", "Crystal Field Theory (CFT)"],
                "topics": ["Ligand: lone pair donor; Lewis base", "Denticity: monodentate, bidentate, polydentate, chelate", "Coordination number: total σ bonds from ligands to metal", "IUPAC: ligands (alphabetical), then metal, then ox. state", "Geometrical isomerism: cis-trans in square planar and octahedral", "Optical isomerism: non-superimposable mirror images (chirality)", "Ionization isomerism: [Co(NH₃)₅Br]SO₄ vs [Co(NH₃)₅SO₄]Br", "Linkage isomerism: ambidentate ligands (SCN⁻, NO₂⁻)", "VBT: hybridization determines geometry (sp = linear, sp³ = tetrahedral, dsp² = square planar, sp³d² = octahedral)", "CFT: crystal field splitting Δ → determines color and magnetism", "Strong field ligands (CO, CN⁻): large Δ, low spin", "Weak field ligands (F⁻, Cl⁻): small Δ, high spin", "Spectrochemical series"],
                "formulas": ["EAN = Z - oxidation state + 2 × CN (18-electron rule for organometallics)", "Color: light absorbed is complementary to color observed"],
                "applications": ["Haemoglobin: Fe²⁺ complex with O₂", "Cisplatin [Pt(NH₃)₂Cl₂]: cancer drug", "EDTA: chelation therapy"],
                "common_mistakes": ["Wrong IUPAC naming order (ligands first alphabetically, then metal)", "Confusing geometrical and optical isomerism conditions", "Wrong hybridization for coordination number"],
                "importance": "12-15% JEE weightage. IUPAC naming, isomerism, and CFT color are annual JEE topics."
            },
            "Haloalkanes & Haloarenes": {
                "difficulty": "Hard",
                "description": "Haloalkanes are alkyl halides (RX) formed by replacing H with a halogen. Haloarenes are aryl halides (ArX). This chapter covers preparation, properties, and reactions of both. SN1 and SN2 are competing nucleophilic substitution mechanisms, with structure determining which predominates. Haloalkanes also serve as synthetic intermediates for making alcohols, ethers, amines, nitriles, and Grignard reagents. Key named reactions (Finkelstein, Swarts, Wurtz-Fittig) are frequently tested.",
                "roadmap": ["Classification & Nomenclature", "Preparation Methods", "SN1 vs SN2 Mechanisms", "Elimination Reactions (E1/E2)", "Reactions & Named Reactions"],
                "topics": ["IUPAC naming of haloalkanes and haloarenes", "Preparation: from alcohols (HX, SOCl₂, PCl₅), from alkenes (HX/X₂), Sandmeyer reaction", "SN2: bimolecular, backside attack, inversion (Walden inversion), 1° halides prefer", "SN1: unimolecular, carbocation intermediate, racemization, 3° halides prefer, polar protic solvent", "E1 (β-elimination): 3° halides + weak nucleophile/base", "E2 (bimolecular): anti-periplanar, strong base, 2°/3° halides", "Reactivity: RI > RBr > RCl > RF", "Grignard reagent: RX + Mg → RMgX (dry ether)", "Finkelstein reaction: RI exchange with NaI in acetone", "Swarts reaction: RF from RCl + AgF or SbF₃", "Wurtz reaction: RX + 2Na → R-R", "Wurtz-Fittig: ArX + R'X + 2Na → Ar-R'", "Haloarenes: C-X bond has partial double bond character → less reactive toward nucleophilic substitution"],
                "formulas": ["SN2: rate = k[RX][Nu⁻] (bimolecular)", "SN1: rate = k[RX] (unimolecular)"],
                "applications": ["Chloroform (CHCl₃) — anesthetic", "DDT — pesticide (banned due to biomagnification)", "Freons (CFCs) — refrigerants (banned due to ozone depletion)"],
                "common_mistakes": ["Wrong mechanism for 3° vs 1° substrates", "Forgetting Grignard reagent requires anhydrous (no water)", "Assuming haloarenes react like haloalkanes in SN reactions"],
                "importance": "12-15% JEE weightage. SN1/SN2 mechanism is JEE staple — appears every year."
            },
            "Alcohols, Phenols & Ethers": {
                "difficulty": "Hard",
                "description": "Alcohols (ROH) and phenols (ArOH) are hydroxyl compounds with contrasting properties — phenols are much more acidic than alcohols due to resonance stabilization of phenoxide ion. Ethers (R-O-R') are relatively unreactive. Lucas test distinguishes 1°, 2°, 3° alcohols; Reimer-Tiemann reaction gives ortho-hydroxybenzaldehyde from phenol.",
                "roadmap": ["Classification & Nomenclature", "Preparation of Alcohols & Phenols", "Physical Properties — H-bonding", "Chemical Reactions", "Phenol Reactions & Acidity"],
                "topics": ["Alcohols: primary (1°), secondary (2°), tertiary (3°)", "Preparation from alkenes (hydration), haloalkanes, carbonyl reduction", "H-bonding: high boiling point compared to alkanes", "Lucas test: ZnCl₂/HCl — 3° (immediate turbidity), 2° (5 min), 1° (no reaction)", "Dehydration: 1° → E2 (180°C), 2°/3° → E1 (low temp)", "Esterification: R-OH + R'COOH → ester + H₂O", "Oxidation: 1° alcohol → aldehyde → carboxylic acid, 2° → ketone, 3° → cleaves (no mild oxidation)", "Victor Meyer test: 1° → red, 2° → blue, 3° → colorless", "Phenol acidity: more acidic than alcohol (pKa ~10) due to resonance", "Electrophilic substitution at ortho and para positions (+M effect of OH)", "Reimer-Tiemann: CHCl₃ + NaOH + Phenol → o-hydroxybenzaldehyde", "Kolbe's reaction: Na-phenoxide + CO₂ → salicylic acid (aspirin precursor)", "Ethers: cleavage by HI/HBr, Williamson synthesis"],
                "formulas": ["Williamson synthesis: R-O⁻ + R'X → R-O-R' (SN2)"],
                "applications": ["Ethanol as biofuel", "Phenol in disinfectants (Dettol)", "Aspirin from salicylic acid (from Kolbe's reaction)"],
                "common_mistakes": ["Forgetting Lucas test only works for tertiary immediately", "Wrong product for Victor Meyer test", "Confusing Reimer-Tiemann with Kolbe reaction"],
                "importance": "15-18% JEE weightage. One of the highest scoring organic chapters."
            },
            "Aldehydes, Ketones & Carboxylic Acids": {
                "difficulty": "Hard",
                "description": "Carbonyl compounds (C=O) are among the most reactive functional groups. Aldehydes (RCHO) are more reactive than ketones (RCOR') toward nucleophilic addition because the absence of an alkyl group makes the carbonyl carbon more electrophilic. Carboxylic acids (RCOOH) are Brønsted acids and undergo esterification, reduction, and halogenation. Named reactions in this chapter are extremely important for JEE.",
                "roadmap": ["Carbonyl Group — Structure & Reactivity", "Nucleophilic Addition Reactions", "Oxidation & Reduction", "Named Reactions (Aldol, Cannizzaro, Clemmensen)", "Carboxylic Acids & Derivatives"],
                "topics": ["Carbonyl group: C=O, sp² hybridized, planar", "Preparation: oxidation of alcohols, ozonolysis, Gattermann, Rosenmund, Stephen", "Nucleophilic addition: HCN, RMgX (Grignard), NaHSO₃, H₂O, alcohols → hemiacetals/acetals", "Reactivity: HCHO > RCHO > RCOR'", "Aldol condensation: α-H + NaOH → β-hydroxy carbonyl compound → enone", "Cannizzaro reaction: HCHO + NaOH (no α-H) → alcohol + acid (disproportionation)", "Clemmensen reduction: RCOR' + Zn/Hg + HCl → RCH₂R' (acidic conditions)", "Wolff-Kishner reduction: RCOR' + H₂NNH₂ + KOH → RCH₂R' (basic conditions)", "Haloform reaction: RCOCH₃ + X₂ + NaOH → RCOONa + CHX₃", "Tollens' test (silver mirror): RCHO + [Ag(NH₃)₂]⁺ → RCOO⁻ + Ag (aldehyde only)", "Fehling's test: RCHO + Cu²⁺ → brick red ppt (aliphatic aldehydes)", "Carboxylic acid: RCOOH, strong H-bonding, high BP", "Esterification: Fischer esterification (acid + alcohol, H⁺ catalyst)", "Hell-Volhard-Zelinsky (HVZ): α-halogenation of carboxylic acid"],
                "formulas": ["Aldol: 2CH₃CHO → CH₃CH(OH)CH₂CHO", "HVZ: RCOOH + X₂/P → RCHXcooH"],
                "applications": ["Formaldehyde (HCHO) in formalin (preservative)", "Acetic acid (CH₃COOH) in vinegar", "Aspirin from salicylic acid esterification"],
                "common_mistakes": ["Cannizzaro reaction requires no α-H — else Aldol occurs", "Tollens' test works for all aldehydes; Fehling's only for aliphatic", "Confusing Clemmensen (acid) with Wolff-Kishner (base)"],
                "importance": "18-20% JEE weightage. Highest scoring organic chapter — named reactions are mandatory."
            },
            "Amines": {
                "difficulty": "Hard",
                "description": "Amines are nitrogen-containing organic bases derived from ammonia. They are classified as primary (1°), secondary (2°), and tertiary (3°) based on the number of alkyl/aryl groups on N. Basicity of aliphatic amines > ammonia > aniline (aromatic amine). Diazonium salts formed from primary aryl amines are key synthetic intermediates for azo dyes, Sandmeyer reactions, and coupling reactions.",
                "roadmap": ["Classification & Nomenclature", "Preparation Methods", "Basicity of Amines", "Reactions", "Diazonium Salts & Coupling"],
                "topics": ["1° (RNH₂), 2° (R₂NH), 3° (R₃N) amines", "Preparation: reduction of nitro compounds, alkylation, Gabriel phthalimide (1° amines only), Hoffmann bromamide (R→R-1 carbon)", "Basicity: 2° aliphatic > 1° aliphatic > 3° aliphatic > NH₃ > aniline", "Aniline: weak base due to lone pair delocalization in ring", "Hinsberg's test: distinguish 1°, 2°, 3° amines using benzene sulfonyl chloride", "Carbylamine test: 1° amines + CHCl₃ + NaOH → isocyanide (foul smell)", "Diazonium salts: ArNH₂ + NaNO₂ + HCl (0-5°C) → ArN₂⁺Cl⁻", "Sandmeyer reaction: ArN₂⁺ + CuCN → ArCN, + CuCl → ArCl, + CuBr → ArBr", "Gattermann reaction: HBF₄ → ArF, HI → ArI", "Coupling reaction: ArN₂⁺ + ArNR₂ → azo dye (orange-red)"],
                "formulas": ["Basicity: 2°(aliphatic) > 1° > 3° > NH₃ > aniline"],
                "applications": ["Aniline — precursor for dyes, drugs, explosives", "Diazonium coupling — azo dyes in textiles", "Lidocaine (local anesthetic) — amide from amine"],
                "common_mistakes": ["Wrong basicity order (it's counterintuitive — 2° > 1° > 3° in gas phase but reversed in water for some)", "Forgetting carbylamine test is only for primary amines", "Diazonium salt preparation requires 0-5°C — higher temp decomposes it"],
                "importance": "12-15% JEE weightage. Named reactions (Hoffmann, Sandmeyer, carbylamine) appear every year."
            },
            "Biomolecules": {
                "difficulty": "Medium",
                "description": "Biomolecules are the large organic molecules that form the basis of all living systems. Carbohydrates provide energy and structural support (cellulose), proteins are functional molecules (enzymes), nucleic acids (DNA, RNA) store and transmit genetic information, and lipids form cell membranes and store energy. This chapter bridges chemistry and biology — important for both boards and JEE.",
                "roadmap": ["Carbohydrates — Classification & Glucose", "Proteins — Amino Acids & Peptides", "Nucleic Acids — DNA & RNA", "Vitamins & Enzymes", "Lipids"],
                "topics": ["Carbohydrates: (CH₂O)n — aldoses and ketoses", "Monosaccharides: glucose (C₆H₁₂O₆), fructose (C₆H₁₂O₆)", "Disaccharides: sucrose (glucose + fructose, non-reducing), maltose/lactose (reducing)", "Polysaccharides: starch (food), cellulose (structural), glycogen (animal storage)", "Open chain structure of glucose: D-glucose, CHO at C-1, OH at C-5 on right", "Haworth structure: pyranose (6-membered) ring form", "Reducing sugars: have free -CHO or -C=O, give Fehling's/Tollens' test", "α-amino acids: H₂N-CHR-COOH, 20 types in proteins", "Peptide bond: -CO-NH- linkage; dipeptide, tripeptide, polypeptide", "Primary structure: sequence of amino acids", "Secondary: α-helix, β-pleated sheet (H-bonds)", "Tertiary: 3D folding; Quaternary: multiple subunits", "Denaturation: loss of structure without breaking peptide bonds", "DNA: double helix (Watson-Crick), A=T (2 H-bonds), G≡C (3 H-bonds)", "RNA: single strand, uracil replaces thymine", "Vitamins: fat-soluble (A, D, E, K) and water-soluble (B complex, C)"],
                "formulas": ["Peptide bond: -CO-NH-", "n amino acids → n-1 peptide bonds + n-1 water molecules"],
                "applications": ["Insulin (protein hormone) for diabetes treatment", "DNA fingerprinting in forensics", "Glycogen as energy reserve in liver"],
                "common_mistakes": ["Sucrose is non-reducing (no free OH on anomeric C)", "Denaturation ≠ hydrolysis of peptide bonds", "Confusing RNA (uracil) with DNA (thymine)"],
                "importance": "8-10% JEE weightage. Glucose structure, DNA base pairing, and reducing sugar tests are frequent."
            },
            "Polymers": {
                "difficulty": "Easy",
                "description": "Polymers are giant molecules (macromolecules) formed by linking many small repeating units called monomers. Addition polymerization involves alkenes (no byproduct), while condensation polymerization forms a small molecule (water, HCl) as byproduct. Natural rubber is a polymer of isoprene; vulcanization adds sulfur crosslinks for strength. Biodegradable polymers (PHBV, nylon-2-nylon-6) are environmentally important.",
                "roadmap": ["Classification of Polymers", "Addition Polymerization (Chain Growth)", "Condensation Polymerization (Step Growth)", "Natural & Synthetic Rubber", "Biodegradable Polymers"],
                "topics": ["Homopolymer (one monomer) vs copolymer (two monomers)", "Linear, branched, cross-linked structures", "Addition polymers: polythene (ethylene), PVC (vinyl chloride), teflon (tetrafluoroethylene), polystyrene, natural rubber", "HDPE (high density, linear) vs LDPE (low density, branched)", "Condensation polymers: nylon-6,6 (hexamethylenediamine + adipic acid), nylon-6 (caprolactam), terylene/dacron (ethylene glycol + terephthalic acid)", "Bakelite: phenol + formaldehyde, thermosetting, cross-linked", "Melamine: melamine + formaldehyde", "Natural rubber: cis-1,4-polyisoprene (2-methyl-1,3-butadiene)n", "Vulcanization: S₈ creates crosslinks, improves elasticity and durability", "Buna-S: butadiene + styrene (styrene-butadiene rubber, SBR)", "Buna-N: butadiene + acrylonitrile (oil-resistant)", "Neoprene: 2-chloro-1,3-butadiene", "Biodegradable: PHBV (3-hydroxybutanoic acid + 3-hydroxypentanoic acid), nylon-2-nylon-6"],
                "formulas": ["Degree of polymerization DP = MW(polymer)/MW(monomer)"],
                "applications": ["LDPE in plastic bags", "Nylon-6,6 in parachutes and stockings", "Bakelite in electrical insulators"],
                "common_mistakes": ["Confusing nylon-6 (one monomer) with nylon-6,6 (two monomers)", "Terylene/dacron is a polyester, not polyamide", "Natural rubber is soft — vulcanization (S) makes it hard"],
                "importance": "6-8% JEE weightage. Monomers of common polymers and Bakelite structure are tested."
            },
            "Chemistry in Everyday Life": {
                "difficulty": "Easy",
                "description": "Chemistry in Everyday Life explores practical applications of chemistry in medicine, food, and hygiene. Medicinal chemistry covers analgesics, antipyretics, antibiotics, antiseptics, disinfectants, antacids, and tranquilizers. Food chemistry covers preservatives, artificial sweeteners, and antioxidants. Cleansing agents cover the mechanism of soaps and detergents. Mostly factual and conceptual — easy scoring chapter for JEE Main.",
                "roadmap": ["Drugs & Medicines Classification", "Analgesics, Antipyretics & Antibiotics", "Antiseptics & Disinfectants", "Food Chemicals", "Cleansing Agents — Soaps & Detergents"],
                "topics": ["Analgesics: non-narcotic (aspirin, paracetamol) and narcotic (morphine, codeine)", "Antipyretics: reduce fever — aspirin, paracetamol", "Antibiotics: kill/inhibit bacteria — penicillin (narrow spectrum), chloramphenicol (broad spectrum)", "Antiseptics: safe for skin (0.2% phenol, iodoform, boric acid)", "Disinfectants: for inanimate objects (1% phenol, bleach)", "Antacids: reduce stomach acidity (baking soda NaHCO₃, milk of magnesia Mg(OH)₂)", "Antihistamines: block histamine action for allergies (diphenhydramine)", "Tranquilizers: treat anxiety and mental disorders (diazepam/Valium, meprobamate)", "Saccharin: artificial sweetener (550× sweeter than sugar), no calories", "Aspartame: sweeter than sugar, breaks down at high temperature", "BHA, BHT: antioxidants in food", "Soap: sodium salt of long chain fatty acid (RCOONa), saponification", "Detergent: sulfonate or sulfate group, works in hard water (no scum)", "Micelle: hydrophobic tail inward, hydrophilic head outward — traps grease"],
                "formulas": ["Saponification: fat + NaOH → soap + glycerol"],
                "applications": ["Aspirin: pain + fever + anti-clotting (prevents heart attacks)", "Penicillin: first antibiotic — from Penicillium mold", "Detergents in hard water (soap forms scum with Ca²⁺/Mg²⁺)"],
                "common_mistakes": ["Antiseptic ≠ disinfectant (antiseptics are safer for living tissue)", "Soap doesn't work in hard water — detergent does", "Saccharin has no food value — aspartame degrades at high temp"],
                "importance": "5-7% JEE weightage. One-liner theory questions — easy scoring in JEE Main."
            },
            "Solutions": {
                "difficulty": "Hard",
                "description": "Solutions explores the properties of liquid mixtures. Concentration expressions (molarity, molality, mole fraction) and Henry's Law form the foundation. Raoult's Law and vapour pressure relate to ideal solutions. Colligative properties — relative lowering of vapour pressure, depression of freezing point (ΔTf), elevation of boiling point (ΔTb), and osmotic pressure (π) — depend only on the number of solute particles (not their nature). Van't Hoff factor (i) accounts for association and dissociation.",
                "roadmap": ["Concentration Expressions", "Henry's Law & Raoult's Law", "Ideal & Non-Ideal Solutions", "Colligative Properties (ΔTf, ΔTb, π)", "Van't Hoff Factor"],
                "topics": ["Molarity M = moles of solute/L of solution", "Molality m = moles of solute/kg of solvent", "Mole fraction X = moles of component/total moles", "Henry's Law: p = KH × x (gas dissolved in liquid)", "Raoult's Law: p = p° × mole fraction of solvent", "Ideal solutions: ΔHmix = 0, ΔVmix = 0 (follow Raoult's over all concentration)", "Positive deviation: A-B weaker than A-A + B-B (e.g., acetone-water)", "Negative deviation: A-B stronger (e.g., acetone-chloroform)", "Azeotropes: constant boiling mixtures at maximum/minimum boiling point", "ΔTf = Kf × m × i (depression of f.p.)", "ΔTb = Kb × m × i (elevation of b.p.)", "π = iMRT (osmotic pressure)", "Van't Hoff factor i: for electrolytes i > 1 (dissociation), for association i < 1", "Osmosis: flow of solvent from low to high concentration through semipermeable membrane", "Reverse osmosis: apply pressure > osmotic pressure → water purification"],
                "formulas": ["ΔTf = Kf × m", "ΔTb = Kb × m", "π = CRT = MRT", "Van't Hoff: i = observed/theoretical colligative property"],
                "applications": ["Antifreeze (ethylene glycol): lowers f.p. of radiator water", "RO water purification plants", "IV saline: isotonic with blood (0.9% NaCl)"],
                "common_mistakes": ["Molarity changes with temperature, molality doesn't", "Forgetting Van't Hoff factor i for electrolytes", "Wrong Raoult's Law for solutions with multiple volatile components"],
                "importance": "12-15% JEE weightage. Colligative property numericals (ΔTf, ΔTb, π) appear every year."
            },
            "Solid State": {
                "difficulty": "Hard",
                "description": "Solid State studies the three-dimensional arrangement of atoms, molecules, or ions in crystalline solids. Crystalline solids have long-range order; amorphous solids do not. Unit cells (cubic: SCC, BCC, FCC) and packing efficiency calculations are heavily tested numerically. Voids (tetrahedral and octahedral), radius ratios, and defects (Schottky and Frenkel) are important conceptual topics.",
                "roadmap": ["Types of Solids — Crystalline vs Amorphous", "Classification of Crystalline Solids", "Unit Cells & Crystal Lattices", "Packing Efficiency & Void Calculation", "Crystal Defects & Properties"],
                "topics": ["Molecular, ionic, covalent, metallic solids", "Crystalline: sharp melting point, anisotropic, long-range order", "Amorphous: no definite melting point, isotropic, short-range order", "7 crystal systems, 14 Bravais lattices", "Simple cubic (SCC): coordination number 6, packing efficiency 52.4%", "BCC: CN = 8, PE = 68%, 2 atoms per unit cell (r = √3a/4)", "FCC/CCP: CN = 12, PE = 74%, 4 atoms per unit cell (r = a/2√2)", "HCP: CN = 12, PE = 74%", "Tetrahedral void: radius ratio 0.225-0.414; Octahedral void: 0.414-0.732", "In FCC: 8 tetrahedral voids, 4 octahedral voids per unit cell", "NaCl structure: FCC with Na in octahedral voids (CN = 6)", "ZnS (zinc blende): FCC with Zn in half tetrahedral voids (CN = 4)", "Schottky defect: equal cation + anion vacancies → decreases density", "Frenkel defect: smaller ion moves to interstitial site → no change in density"],
                "formulas": ["Density = Z × M / (a³ × NA)", "Packing efficiency (FCC) = 74%", "Packing efficiency (BCC) = 68%"],
                "applications": ["Diamond (covalent crystal) hardest solid", "NaCl crystal structure electronics", "Semiconductor doping (crystal defects)"],
                "common_mistakes": ["Counting atoms in unit cell (corner = 1/8, face = 1/2, edge = 1/4, body = 1)", "Schottky vs Frenkel: Schottky = vacancy, Frenkel = interstitial", "Wrong coordination number for different structures"],
                "importance": "10-12% JEE weightage. Unit cell density formula and packing efficiency are annual JEE numericals."
            }
        },

        "Maths": {
            "Calculus": {
                "difficulty": "Hard",
                "roadmap": ["Continuity & Differentiability", "Differentiation Applications", "Integration Techniques", "Definite Integrals", "Area Under Curves"],
                "topics": ["Derivative Rules", "Rolle's & Lagrange's Theorem", "Maxima & Minima", "Integration by Parts", "Definite Integral Properties"],
                "formulas": ["d/dx(xⁿ) = nxⁿ⁻¹", "∫x^n dx = x^(n+1)/(n+1)", "∫udv = uv - ∫vdu", "dy/dx = (dy/dt)/(dx/dt)"],
                "importance": "Highest 25-30% board weightage. Core JEE Advanced topic."
            },
            "Vectors": {
                "difficulty": "Medium",
                "roadmap": ["Vector Algebra Basics", "Dot Product", "Cross Product", "Scalar Triple Product", "Applications"],
                "topics": ["Direction Cosines & Ratios", "Position Vector", "Projections", "Area of Parallelogram", "Scalar & Vector Triple Products"],
                "formulas": ["a·b = |a||b|cosθ", "a×b = |a||b|sinθ n̂", "[a b c] = a·(b×c)"],
                "importance": "10-12% board weightage. 3D geometry foundation for JEE."
            },
            "3D Geometry": {
                "difficulty": "Hard",
                "roadmap": ["Direction Cosines", "Equation of Line", "Equation of Plane", "Angle Between Lines/Planes", "Distance Formulas"],
                "topics": ["Cartesian & Vector Form", "Skew Lines", "Coplanar Lines", "Perpendicular Distance", "Shortest Distance"],
                "formulas": ["(x-x₁)/a = (y-y₁)/b = (z-z₁)/c", "ax + by + cz + d = 0", "d = |ax₁+by₁+cz₁+d|/√(a²+b²+c²)"],
                "importance": "15% board weightage. 3D visualization needed for JEE."
            },
            "Matrices & Determinants": {
                "difficulty": "Hard",
                "roadmap": ["Matrix Operations", "Types of Matrices", "Determinant Properties", "Adjoint & Inverse", "System of Equations"],
                "topics": ["Addition, Multiplication", "Transpose, Symmetric", "Cofactor Expansion", "Cramer's Rule", "Rank of Matrix"],
                "formulas": ["|AB| = |A||B|", "A⁻¹ = adj(A)/|A|", "A(adj A) = |A|I"],
                "importance": "12-15% board weightage. Matrix problems frequent in JEE."
            },
            "Probability": {
                "difficulty": "Hard",
                "roadmap": ["Random Experiments", "Conditional Probability", "Bayes' Theorem", "Random Variables", "Probability Distributions"],
                "topics": ["Addition & Multiplication Theorems", "Independent Events", "Binomial Distribution", "Mean & Variance"],
                "formulas": ["P(A∪B) = P(A) + P(B) - P(A∩B)", "P(A|B) = P(A∩B)/P(B)", "P(E) = ΣP(Ei)P(E|Ei)"],
                "importance": "10-12% board weightage. Combinatorics integration in JEE."
            },
            "Differential Equations": {
                "difficulty": "Hard",
                "roadmap": ["Order & Degree", "Variables Separable", "Homogeneous Equations", "Linear Equations", "Applications"],
                "topics": ["First Order DE", "dy/dx + Py = Q", "Integrating Factor", "Formation of DE", "Growth & Decay"],
                "formulas": ["IF = e^∫Pdx", "y·IF = ∫Q·IF dx"],
                "importance": "8-10% board weightage. Application-based JEE problems."
            },
            "Relations & Functions": {
                "difficulty": "Medium",
                "roadmap": ["Types of Relations", "Types of Functions", "Composition", "Inverse Functions", "Binary Operations"],
                "topics": ["Reflexive, Symmetric, Transitive", "One-One, Onto, Bijective", "fog & gof", "Invertible Functions"],
                "formulas": ["(fog)(x) = f(g(x))", "f⁻¹(f(x)) = x"],
                "importance": "8% board weightage. Conceptual clarity for JEE."
            },
            "Linear Programming": {
                "difficulty": "Easy",
                "roadmap": ["Problem Formulation", "Graphical Solution", "Feasible Region", "Optimization", "Corner Point Method"],
                "topics": ["Objective Function", "Constraints", "Bounded & Unbounded Regions", "Maximum & Minimum Values"],
                "formulas": ["Optimize Z = ax + by subject to constraints"],
                "importance": "6-8% board weightage. Easy graphical marks."
            },
            "Applications of Integrals": {
                "difficulty": "Medium",
                "roadmap": ["Area Between Curves", "Area of Circle/Ellipse", "Volume of Solids"],
                "topics": ["Area under y = f(x)", "Area between two curves", "Volume by rotation"],
                "formulas": ["Area = ∫[a to b] f(x)dx", "Area = ∫[a to b] [f(x) - g(x)]dx"],
                "importance": "8-10% board weightage. Visualization needed for JEE."
            }
        }
    }
}

extra_data = {
    "9": {
        "Physics": {
            "Motion": {
                "difficulty": "Easy",
                "description": "Motion is the change in position of an object with respect to time and a reference point. This chapter introduces students to the concepts of distance, displacement, speed, velocity, and acceleration. Through graphs and equations, students learn to analyze the motion of objects in a straight line and develop problem-solving skills that form the base for all of mechanics.",
                "roadmap": ["Distance vs Displacement", "Speed vs Velocity", "Acceleration", "Equations of Motion", "Distance-Time & Velocity-Time Graphs"],
                "topics": ["Uniform & Non-Uniform Motion", "Scalar & Vector Quantities", "Average Speed & Velocity", "Uniform Acceleration", "Graphical Analysis of Motion"],
                "formulas": ["v = u + at", "s = ut + ½at²", "v² = u² + 2as", "s = (u+v)/2 × t"],
                "importance": "Foundation for kinematics in Class 11. 10-12% board weightage. Easy scoring chapter."
            },
            "Force & Laws of Motion": {
                "difficulty": "Medium",
                "description": "Force & Laws of Motion introduces Newton's three laws that explain why objects move, stop, or change direction. It builds intuition about inertia, momentum, and action-reaction forces through real-life examples such as seat belts, rockets, and sports. This chapter is critical for understanding all of classical mechanics.",
                "roadmap": ["Force & Its Types", "Newton's First Law - Inertia", "Newton's Second Law - F=ma", "Newton's Third Law - Action-Reaction", "Conservation of Momentum"],
                "topics": ["Balanced & Unbalanced Forces", "Types of Inertia", "Linear Momentum", "Impulse", "Rocket Propulsion"],
                "formulas": ["F = ma", "p = mv", "F = Δp/Δt", "m₁u₁ + m₂u₂ = m₁v₁ + m₂v₂"],
                "importance": "Core mechanics concept. 15% board weightage. Base for Class 11 Laws of Motion."
            },
            "Gravitation": {
                "difficulty": "Easy",
                "description": "Gravitation covers the universal law of attraction between masses and explains why objects fall towards Earth. It introduces concepts of free fall, weight vs mass, and the universal gravitational constant G. Students also learn about satellite motion and the conditions required for weightlessness.",
                "roadmap": ["Universal Law of Gravitation", "Free Fall & g", "Weight vs Mass", "Kepler's Laws Introduction", "Satellites & Weightlessness"],
                "topics": ["Gravitational Force", "Factors Affecting g", "Equations for Free Fall", "Escape Velocity Basics", "Moon's Motion"],
                "formulas": ["F = Gm₁m₂/r²", "g = GM/R²", "W = mg", "v = u + gt", "h = ut + ½gt²"],
                "importance": "8-10% board weightage. Conceptual understanding needed for Class 11 Gravitation."
            },
            "Work & Energy": {
                "difficulty": "Easy",
                "description": "Work and Energy explains how energy is transferred when a force causes displacement. It introduces kinetic and potential energy, and establishes the work-energy theorem and law of conservation of energy. Power is defined as the rate of doing work, and these concepts are applied to machines and everyday devices.",
                "roadmap": ["Definition of Work", "Kinetic Energy", "Potential Energy", "Conservation of Energy", "Power & Commercial Unit"],
                "topics": ["Conditions for Work", "Work Done by Gravity", "Elastic & Gravitational PE", "Energy Transformation", "Commercial Unit of Energy (kWh)"],
                "formulas": ["W = F × s × cosθ", "KE = ½mv²", "PE = mgh", "Power = W/t", "1 kWh = 3.6 × 10⁶ J"],
                "importance": "10-12% board weightage. Easy numericals. Foundation for Class 11 Work Energy Power."
            },
            "Sound": {
                "difficulty": "Medium",
                "description": "Sound explores the nature of sound waves as mechanical longitudinal waves that require a medium for propagation. The chapter covers the characteristics of sound such as pitch, loudness, and quality, as well as important phenomena like reflection (echo), reverberation, and the SONAR technique used in underwater exploration.",
                "roadmap": ["Nature of Sound Waves", "Propagation of Sound", "Characteristics of Sound", "Reflection of Sound", "Human Ear & SONAR"],
                "topics": ["Compression & Rarefaction", "Speed of Sound in Different Media", "Frequency, Wavelength, Amplitude", "Echo & Reverberation", "Ultrasound Applications"],
                "formulas": ["v = fλ", "Speed of sound in air ≈ 344 m/s at 22°C", "Minimum distance for echo = 17.2 m"],
                "importance": "10% board weightage. Connects to Waves chapter in Class 11."
            }
        },
        "Chemistry": {
            "Matter in Our Surroundings": {
                "difficulty": "Easy",
                "description": "Matter in Our Surroundings introduces students to the physical nature of matter — anything that has mass and occupies space. It discusses the three states of matter (solid, liquid, gas) and explains changes of state including evaporation, boiling, melting, and sublimation based on particle behavior and intermolecular forces.",
                "roadmap": ["States of Matter", "Characteristics of Particles", "Change of State", "Evaporation", "Latent Heat"],
                "topics": ["Solid, Liquid, Gas Properties", "Interconversion of States", "Boiling vs Evaporation", "Sublimation Examples", "Effect of Temperature & Pressure"],
                "formulas": ["No major formulas - Conceptual"],
                "importance": "8% board weightage. Easy scoring introductory chapter."
            },
            "Atoms & Molecules": {
                "difficulty": "Medium",
                "description": "Atoms & Molecules establishes the particulate nature of matter. It covers Dalton's atomic theory, laws of chemical combination, concept of mole, and calculation of molar mass and number of particles. This chapter builds the quantitative foundation of chemistry.",
                "roadmap": ["Dalton's Atomic Theory", "Laws of Chemical Combination", "Atoms & Symbols", "Molecules & Formulae", "Mole Concept"],
                "topics": ["Law of Conservation of Mass", "Law of Definite Proportions", "Atomic Mass & AMU", "Molecular Mass", "Avogadro's Number"],
                "formulas": ["Molar mass = mass of 6.022×10²³ particles", "No. of moles = given mass / molar mass", "No. of particles = moles × 6.022×10²³"],
                "importance": "12-15% board weightage. Foundation for stoichiometry."
            },
            "Structure of Atom": {
                "difficulty": "Medium",
                "description": "Structure of Atom traces the historical development of atomic models — from Thomson's plum pudding model to Rutherford's nuclear model to Bohr's model. Students learn about subatomic particles, electron distribution in shells, valency, and atomic number vs mass number.",
                "roadmap": ["Thomson's Model", "Rutherford's Experiment", "Bohr's Model", "Electronic Configuration", "Valency & Atomic Number"],
                "topics": ["Protons, Neutrons, Electrons", "Atomic Number & Mass Number", "Isotopes & Isobars", "Bohr's Shells (K,L,M,N)", "2n² Rule"],
                "formulas": ["Mass Number = Protons + Neutrons", "Electrons in shell = 2n²", "Valency from configuration"],
                "importance": "12% board weightage. Foundation for Class 11 Atomic Structure."
            },
            "Chemical Reactions & Equations": {
                "difficulty": "Easy",
                "description": "Chemical Reactions & Equations teaches students to write and balance chemical equations using the law of conservation of mass. It classifies reactions into combination, decomposition, displacement, double displacement, oxidation, and reduction — and explores everyday applications like rusting and rancidity.",
                "roadmap": ["Writing Chemical Equations", "Balancing Equations", "Types of Reactions", "Oxidation & Reduction", "Corrosion & Rancidity"],
                "topics": ["Combination Reactions", "Decomposition Reactions", "Displacement Reactions", "Redox Reactions", "Effects of Corrosion"],
                "formulas": ["Law of Conservation of Mass", "Balancing by Hit & Trial"],
                "importance": "10-12% board weightage. Practical application-based chapter."
            },
            "Acids, Bases & Salts": {
                "difficulty": "Medium",
                "description": "Acids, Bases & Salts introduces the chemical nature of common substances. It explains acid-base theories, pH scale, neutralization reactions, and preparation of salts. Everyday examples like baking soda, washing soda, bleaching powder, and plaster of paris make this chapter highly practical and interesting.",
                "roadmap": ["Acids & Their Properties", "Bases & Their Properties", "Neutralization", "pH Scale", "Salts & Their Types"],
                "topics": ["Strong & Weak Acids/Bases", "Universal Indicator", "Antacids & Baking Soda", "Washing Soda vs Baking Soda", "Bleaching Powder"],
                "formulas": ["pH = -log[H⁺]", "Acid + Base → Salt + Water", "pH < 7 Acidic, pH > 7 Basic"],
                "importance": "15% board weightage. High-scoring practical chapter."
            }
        },
        "Maths": {
            "Number Systems": {
                "difficulty": "Easy",
                "description": "Number Systems explores the hierarchy of numbers — natural numbers, whole numbers, integers, rational and irrational numbers — and their representation on the number line. Students learn to represent surds on number lines, rationalize denominators, and apply exponent rules. This chapter strengthens numerical foundations for all of mathematics.",
                "roadmap": ["Natural, Whole, Integer, Rational Numbers", "Irrational Numbers", "Real Numbers on Number Line", "Surds Representation", "Laws of Exponents"],
                "topics": ["Rational Number as p/q", "Non-terminating Recurring Decimals", "Representation of √2 on Number Line", "Rationalisation", "Laws of Radicals"],
                "formulas": ["(√a + √b)(√a - √b) = a - b", "(a + b)² = a² + 2ab + b²", "aᵐ × aⁿ = aᵐ⁺ⁿ"],
                "importance": "8-10% board weightage. Strong base for algebra and trigonometry."
            },
            "Polynomials": {
                "difficulty": "Medium",
                "description": "Polynomials introduces algebraic expressions with one or more variables and covers types, degrees, zeroes, and operations on polynomials. The Remainder and Factor theorems simplify solving polynomial equations. Algebraic identities are tools that appear across the entire mathematics curriculum.",
                "roadmap": ["Types & Degrees of Polynomials", "Zeroes of Polynomial", "Remainder Theorem", "Factor Theorem", "Algebraic Identities"],
                "topics": ["Monomial, Binomial, Trinomial", "Division Algorithm", "Factorisation using Identities", "Relationship between Zeroes & Coefficients"],
                "formulas": ["(a+b)³ = a³+3a²b+3ab²+b³", "(a-b)³ = a³-3a²b+3ab²-b³", "a³+b³ = (a+b)(a²-ab+b²)"],
                "importance": "12% board weightage. Algebra skills crucial for JEE prep."
            },
            "Triangles": {
                "difficulty": "Medium",
                "description": "Triangles covers congruence and similarity criteria for triangles, including SAS, ASA, SSS, RHS, and AA. The chapter proves the Pythagoras theorem and its converse, and introduces the Basic Proportionality Theorem. These concepts are the foundation of Euclidean geometry.",
                "roadmap": ["Congruence Criteria", "Similarity Criteria", "Pythagoras Theorem", "Basic Proportionality Theorem", "Areas of Similar Triangles"],
                "topics": ["SAS, ASA, SSS, RHS", "AA, SAS, SSS Similarity", "Proof of Pythagoras", "Ratio of Areas of Similar Triangles"],
                "formulas": ["a² + b² = c² (Pythagoras)", "Ratio of areas = (ratio of sides)²", "BPT: DE ∥ BC → AD/DB = AE/EC"],
                "importance": "15% board weightage. Geometry base for higher classes."
            },
            "Coordinate Geometry": {
                "difficulty": "Easy",
                "description": "Coordinate Geometry introduces the Cartesian plane and teaches students to locate points, find distances between points, and calculate the coordinates of midpoints and section points. It bridges algebra and geometry by showing how geometric relationships can be expressed as algebraic equations.",
                "roadmap": ["Cartesian Plane", "Distance Formula", "Section Formula", "Midpoint Formula", "Area of Triangle"],
                "topics": ["Plotting Points in Four Quadrants", "Distance Between Two Points", "Internal & External Division", "Collinearity Check", "Area of Triangle on Coordinate Plane"],
                "formulas": ["Distance = √[(x₂-x₁)² + (y₂-y₁)²]", "Midpoint = ((x₁+x₂)/2, (y₁+y₂)/2)", "Area = ½|x₁(y₂-y₃)+x₂(y₃-y₁)+x₃(y₁-y₂)|"],
                "importance": "10-12% board weightage. Connects to straight lines in Class 11."
            },
            "Statistics": {
                "difficulty": "Easy",
                "description": "Statistics introduces the collection, organization, presentation, and analysis of numerical data. Students learn to calculate mean, median, and mode for ungrouped and grouped data, and represent data using bar graphs, histograms, and frequency polygons. These tools are essential for data literacy in everyday life.",
                "roadmap": ["Data Collection & Organisation", "Bar Graphs & Histograms", "Frequency Polygons", "Mean, Median, Mode", "Central Tendency Applications"],
                "topics": ["Primary vs Secondary Data", "Frequency Distribution Table", "Class Intervals", "Arithmetic Mean by Direct Method", "Mode of Grouped Data"],
                "formulas": ["Mean = Σfx/Σf", "Median for odd n = ((n+1)/2)th term", "Mode = most frequent value"],
                "importance": "10% board weightage. Easy scoring and practically relevant."
            }
        }
    },
    "10": {
        "Physics": {
            "Electricity": {
                "difficulty": "Medium",
                "description": "Electricity introduces electric charge, current, potential difference, resistance, and Ohm's law. Students learn about factors affecting resistance, series and parallel circuits, and heating effects of current. This chapter forms the practical core of Class 10 Physics and is heavily tested in board exams with numerical problems.",
                "roadmap": ["Electric Charge & Current", "Potential Difference & Ohm's Law", "Resistance & Factors", "Series & Parallel Circuits", "Heating Effect & Power"],
                "topics": ["Coulomb, Ampere, Volt, Ohm", "V-I Characteristics", "Resistivity", "Kirchhoff's Rules (Basic)", "Joule's Law of Heating"],
                "formulas": ["V = IR", "R = ρl/A", "Rs = R₁ + R₂", "1/Rp = 1/R₁ + 1/R₂", "P = VI = I²R = V²/R", "H = I²Rt"],
                "importance": "20% board weightage. Maximum numericals. Critical for Class 12 Current Electricity."
            },
            "Magnetic Effects of Electric Current": {
                "difficulty": "Medium",
                "description": "Magnetic Effects of Electric Current covers how electric currents produce magnetic fields and how magnets interact with current-carrying conductors. It explains the working of electric motors, electromagnetic induction, and electric generators — the technology behind our modern energy infrastructure.",
                "roadmap": ["Magnetic Field of a Current", "Electromagnets", "Force on Current in Magnetic Field", "Electric Motor", "Electromagnetic Induction & Generator"],
                "topics": ["Right-Hand Thumb Rule", "Magnetic Field of Solenoid", "Fleming's Left-Hand Rule", "DC Motor Working", "AC vs DC Generator"],
                "formulas": ["F = BIlsinθ", "Electromagnetic Induction - Faraday's Law (qualitative)"],
                "importance": "15-18% board weightage. Conceptual + diagram-heavy chapter."
            },
            "Light - Reflection & Refraction": {
                "difficulty": "Medium",
                "description": "Light - Reflection & Refraction covers the laws of reflection and refraction and their applications in mirrors and lenses. Sign conventions, ray diagrams, and numerical problems on mirror and lens formula, magnification, and power of lenses are key skills students must master in this chapter.",
                "roadmap": ["Laws of Reflection", "Spherical Mirrors", "Refraction Laws", "Spherical Lenses", "Power of Lens"],
                "topics": ["Concave & Convex Mirrors", "Mirror Formula & Magnification", "Snell's Law", "Concave & Convex Lens", "Human Eye & Defects"],
                "formulas": ["1/f = 1/v + 1/u", "m = -v/u = h'/h", "n = sin i / sin r", "P = 1/f (in metres)", "P_total = P₁ + P₂"],
                "importance": "20% board weightage. Maximum diagram marks. Foundation for Class 12 Optics."
            },
            "Sources of Energy": {
                "difficulty": "Easy",
                "description": "Sources of Energy surveys conventional and non-conventional energy sources. It discusses fossil fuels, thermal and hydroelectric power, nuclear energy, solar energy, wind energy, and biomass. The chapter emphasizes the environmental impacts of energy use and the need for sustainable, renewable alternatives.",
                "roadmap": ["Conventional Energy Sources", "Thermal Power", "Hydroelectric Power", "Non-Conventional Sources", "Nuclear Energy"],
                "topics": ["Fossil Fuels & Limitations", "Biogas Plants", "Solar Cookers & Cells", "Wind Energy", "Nuclear Fission & Fusion Basics"],
                "formulas": ["Efficiency = Useful Output / Total Input × 100"],
                "importance": "8-10% board weightage. Theory-based easy marks."
            }
        },
        "Chemistry": {
            "Chemical Reactions & Equations": {
                "difficulty": "Easy",
                "description": "Building on Class 9 basics, this chapter focuses on writing and balancing chemical equations with state symbols and conditions. It classifies reactions as combination, decomposition, displacement, double displacement, and redox — with practical applications including rusting, rancidity, and precipitation reactions.",
                "roadmap": ["Balanced Equations with State Symbols", "Types of Chemical Reactions", "Oxidation & Reduction (Redox)", "Corrosion & Rancidity", "Exothermic & Endothermic Reactions"],
                "topics": ["Balancing Methods", "Combination & Decomposition", "Single & Double Displacement", "Thermal Decomposition", "Precipitation Reactions"],
                "formulas": ["2H₂ + O₂ → 2H₂O", "CaCO₃ → CaO + CO₂", "Fe + CuSO₄ → FeSO₄ + Cu"],
                "importance": "12-15% board weightage. Foundation of inorganic chemistry."
            },
            "Acids, Bases & Salts": {
                "difficulty": "Medium",
                "description": "This chapter deepens the understanding of acids and bases through Arrhenius and Brønsted-Lowry definitions. It covers the pH scale in detail, common salts like baking soda, washing soda, bleaching powder, and plaster of paris — their preparation, properties, and uses.",
                "roadmap": ["Acids Properties & Reactions", "Bases Properties & Reactions", "Salts - Neutral, Acidic, Basic", "pH Calculations", "Important Salts & Uses"],
                "topics": ["Dilution of Acids", "Metal Oxides vs Non-metal Oxides", "Salt Hydrolysis", "Sodium Hydroxide Production", "Water of Crystallisation"],
                "formulas": ["pH = -log[H⁺]", "pH + pOH = 14", "NaCl + H₂O + CO₂ + NH₃ → NaHCO₃ + NH₄Cl (Solvay Process)"],
                "importance": "15-18% board weightage. High scoring chapter with easy numericals."
            },
            "Metals & Non-Metals": {
                "difficulty": "Medium",
                "description": "Metals & Non-Metals compares physical and chemical properties of metals and non-metals. It explains the reactivity series, extraction of metals (metallurgy), formation of ionic bonds, and corrosion prevention. Understanding this chapter is crucial for inorganic chemistry at higher levels.",
                "roadmap": ["Physical Properties Comparison", "Chemical Properties of Metals", "Reactivity Series", "Extraction of Metals", "Corrosion & Prevention"],
                "topics": ["Malleability, Ductility, Conductivity", "Reaction with Water, Acid, Oxygen", "Highly Reactive to Noble Metals", "Reduction, Refining, Ore Concentration", "Galvanisation & Alloying"],
                "formulas": ["2Na + 2H₂O → 2NaOH + H₂", "Fe₂O₃ + 3CO → 2Fe + 3CO₂", "2Al₂O₃ → 4Al + 3O₂ (Electrolysis)"],
                "importance": "15% board weightage. Inorganic base for JEE preparation."
            },
            "Carbon & Its Compounds": {
                "difficulty": "Medium",
                "description": "Carbon & Its Compounds introduces organic chemistry through the unique bonding properties of carbon — its ability to form chains, rings, and multiple bonds. Students learn about homologous series, functional groups, IUPAC naming, and important organic compounds like alkanes, alkenes, ethanol, and ethanoic acid.",
                "roadmap": ["Covalent Bonding in Carbon", "Allotropes of Carbon", "Hydrocarbons", "Functional Groups & Homologous Series", "Ethanol & Ethanoic Acid"],
                "topics": ["Tetravalency & Catenation", "Diamond, Graphite, Fullerene", "Saturated & Unsaturated Hydrocarbons", "Nomenclature", "Chemical Properties of Ethanol & Ethanoic Acid"],
                "formulas": ["CₙH₂ₙ₊₂ (Alkanes)", "CₙH₂ₙ (Alkenes)", "CₙH₂ₙ₋₂ (Alkynes)", "CH₃COOH (Acetic Acid)"],
                "importance": "18-20% board weightage. Maximum organic marks in Class 10."
            },
            "Periodic Classification of Elements": {
                "difficulty": "Easy",
                "description": "Periodic Classification traces the history of element organization — from Döbereiner's triads to Newlands' octaves to Mendeleev's periodic table and the modern long form. Students understand how the periodic table predicts element properties and how electronic configuration determines an element's position.",
                "roadmap": ["Döbereiner's Triads", "Newlands' Law of Octaves", "Mendeleev's Periodic Table", "Modern Periodic Law", "Periodic Trends"],
                "topics": ["Merits & Limitations of Each Table", "Groups & Periods", "Atomic Size, Metallic Character", "Valency Trends", "Position & Electronic Configuration"],
                "formulas": ["Group = number of valence electrons", "Period = number of shells"],
                "importance": "10% board weightage. Conceptual clarity for Class 11 Periodic Properties."
            }
        },
        "Maths": {
            "Real Numbers": {
                "difficulty": "Easy",
                "description": "Real Numbers starts with Euclid's Division Lemma — a powerful tool to find the HCF of two numbers through repeated division. The Fundamental Theorem of Arithmetic establishes that every composite number has a unique prime factorization. The chapter also proves the irrationality of numbers like √2 and √3 and explores decimal expansions.",
                "roadmap": ["Euclid's Division Lemma", "Fundamental Theorem of Arithmetic", "Prime Factorisation", "Irrational Numbers Proof", "Decimal Expansions"],
                "topics": ["HCF by Euclid's Algorithm", "LCM × HCF = Product of Numbers", "Proof that √2 is Irrational", "Terminating vs Non-terminating Decimals"],
                "formulas": ["a = bq + r (Euclid's Lemma)", "HCF × LCM = a × b", "p/q is terminating if q = 2ᵐ × 5ⁿ"],
                "importance": "6-8% board weightage. Conceptual base for number theory."
            },
            "Polynomials": {
                "difficulty": "Medium",
                "description": "Polynomials in Class 10 focuses on the relationship between the zeroes of quadratic and cubic polynomials and their coefficients. Students verify zeroes, find polynomials given zeroes, and divide polynomials using the division algorithm. These skills underpin factorization across all of algebra.",
                "roadmap": ["Geometrical Meaning of Zeroes", "Relationship Between Zeroes & Coefficients", "Division Algorithm", "Finding Polynomial from Zeroes"],
                "topics": ["Linear, Quadratic, Cubic Polynomials", "Sum & Product of Zeroes", "Verification of Zeroes", "Long Division of Polynomials"],
                "formulas": ["Sum of zeroes α+β = -b/a", "Product of zeroes αβ = c/a", "For cubic: α+β+γ = -b/a, αβ+βγ+γα = c/a, αβγ = -d/a"],
                "importance": "10-12% board weightage. Algebraic skills for Class 11."
            },
            "Quadratic Equations": {
                "difficulty": "Medium",
                "description": "Quadratic Equations covers solving second-degree equations by factorization, completing the square, and the quadratic formula. The discriminant determines the nature of roots. Word problems involving area, speed, time, and age make this chapter highly practical and frequently tested.",
                "roadmap": ["Standard Form", "Solution by Factorisation", "Completing the Square", "Quadratic Formula", "Nature of Roots via Discriminant"],
                "topics": ["Roots by Factorisation", "Quadratic Formula Derivation", "Discriminant D = b² - 4ac", "Real vs Complex Roots", "Word Problems"],
                "formulas": ["ax² + bx + c = 0", "x = [-b ± √(b²-4ac)] / 2a", "D = b² - 4ac", "Sum = -b/a, Product = c/a"],
                "importance": "12-15% board weightage. High-scoring chapter with consistent patterns."
            },
            "Triangles": {
                "difficulty": "Medium",
                "description": "Triangles in Class 10 extends Class 9 concepts to similarity with the Basic Proportionality Theorem and its converse. Students prove and apply similarity criteria (AA, SSS, SAS), establish the relationship between areas of similar triangles, and extend Pythagoras theorem. This chapter has the highest theorem-proof density in Class 10 Maths.",
                "roadmap": ["Basic Proportionality Theorem", "Criteria for Similarity", "Areas of Similar Triangles", "Pythagoras Theorem", "Converse of Pythagoras"],
                "topics": ["Thales' Theorem", "AA Similarity", "Ratio of Areas = Square of Ratio of Sides", "Proof of Pythagoras", "Acute & Obtuse Triangles"],
                "formulas": ["AD/DB = AE/EC (BPT)", "Area ratio = (side ratio)²", "a² + b² = c²", "a² = b² + c² - 2bc cosA"],
                "importance": "15-20% board weightage. Highest proof marks. Very important."
            },
            "Statistics & Probability": {
                "difficulty": "Medium",
                "description": "Statistics covers mean, median, and mode for grouped data using direct, assumed mean, and step deviation methods. Cumulative frequency curves (ogives) are drawn to estimate median graphically. Probability introduces the classical definition and calculates theoretical probability for simple experiments involving coins, dice, and cards.",
                "roadmap": ["Mean of Grouped Data", "Mode of Grouped Data", "Median of Grouped Data", "Ogive & Cumulative Frequency", "Probability - Classical Definition"],
                "topics": ["Direct Method, Assumed Mean, Step Deviation", "Modal Class", "Less Than & More Than Ogive", "Probability of Simple Events", "Complementary Events"],
                "formulas": ["Mean = a + (Σfd/Σf)×h", "Mode = l + ((f₁-f₀)/(2f₁-f₀-f₂))×h", "Median = l + ((n/2 - cf)/f)×h", "P(E) = n(E)/n(S)"],
                "importance": "15% board weightage. Practical and scoring chapter."
            }
        }
    }
}
data.update(extra_data)

data["11"]["Maths"]["Permutations & Combinations"]["description"] = "Permutations & Combinations is the mathematics of counting and arranging. It introduces the fundamental principles of counting — the multiplication and addition rules — and builds up to permutations (ordered arrangements) and combinations (unordered selections). This chapter is crucial for probability and appears frequently in JEE with complex constraint-based problems."
data["11"]["Maths"]["Binomial Theorem"]["description"] = "Binomial Theorem provides a formula to expand expressions of the form (a+b)ⁿ without multiplying out term by term. It introduces binomial coefficients, Pascal's triangle, and methods to find specific terms in the expansion. Applications range from approximations in science to algebraic proofs in higher mathematics."
data["11"]["Maths"]["Trigonometry"]["description"] = "Trigonometry studies the relationships between the sides and angles of triangles. It extends beyond triangles to define functions — sine, cosine, tangent and their reciprocals — across all angles. Mastery of identities, inverse functions, and equations is vital for calculus, coordinate geometry, and physics problems throughout Class 11 and JEE."
data["11"]["Maths"]["Sets & Relations"]["description"] = "Sets & Relations is the language of modern mathematics. Set theory defines collections of distinct objects and their operations (union, intersection, complement). Relations connect elements of different sets, and functions are special relations. This chapter develops logical and abstract thinking needed for all of mathematics."
data["11"]["Maths"]["Sequences & Series"]["description"] = "Sequences & Series studies ordered lists of numbers generated by specific patterns — arithmetic progressions (AP), geometric progressions (GP), and harmonic progressions (HP). It develops tools to find any term, calculate sums, and evaluate infinite geometric series — concepts that appear in calculus through series expansions."
data["11"]["Maths"]["Straight Lines"]["description"] = "Straight Lines extends coordinate geometry to analyze lines algebraically. Students learn to find equations of lines in various forms, calculate angles between lines, determine distances from points to lines, and understand the relationship between parallel and perpendicular lines. This chapter is the gateway to conic sections and 3D geometry."
data["11"]["Maths"]["Conic Sections"]["description"] = "Conic Sections studies the curves formed by intersecting a cone with a plane — circles, parabolas, ellipses, and hyperbolas. Each curve has standard equations, key geometric properties (focus, directrix, eccentricity), and applications in physics (satellite orbits, projectile paths) and engineering (antenna design)."
data["11"]["Maths"]["Limits & Derivatives"]["description"] = "Limits & Derivatives introduces differential calculus — one of the most powerful tools in mathematics. Limits formalize the concept of approaching a value, while derivatives measure instantaneous rates of change. This chapter establishes the foundation for all of Class 12 calculus, including differentiation, integration, and their real-world applications."
data["12"]["Chemistry"]["Polymers"]["description"] = "Polymers are giant molecules formed by the repetitive linking of small molecules called monomers. This chapter classifies polymers by structure, source, and mode of formation — covering addition and condensation polymerization. Students study natural rubber, synthetic polymers like nylon and polyester, and biodegradable alternatives, exploring how molecular structure determines material properties."
data["12"]["Chemistry"]["Biomolecules"]["description"] = "Biomolecules studies the chemistry of life — the large organic molecules that make up living organisms. Carbohydrates provide energy, proteins carry out cellular functions, nucleic acids store genetic information, and lipids form cell membranes. Understanding their structure, classification, and reactions is essential for biochemistry, medicine, and biotechnology."
data["12"]["Chemistry"]["Chemistry in Everyday Life"]["description"] = "Chemistry in Everyday Life connects abstract chemistry to practical applications. It covers medicinal chemistry (analgesics, antibiotics, antiseptics, tranquilizers, antacids), food chemistry (preservatives, sweeteners, antioxidants), and cleansing agents (soaps, detergents). This chapter helps students appreciate how chemistry shapes health, food safety, and hygiene."

@app.route('/')
def home():
    return send_from_directory('.', 'index.html')

@app.route('/api/class/<class_id>/subject/<subject>')
def get_subject_data(class_id, subject):
    return jsonify(data.get(class_id, {}).get(subject, {}))

@app.route('/api/class/<class_id>/subject/<subject>/topic/<topic>')
def get_topic_data(class_id, subject, topic):
    subject_data = data.get(class_id, {}).get(subject, {})
    for key in subject_data:
        if key.lower().replace(" ", "") == topic.lower().replace(" ", ""):
            return jsonify(subject_data[key])
    return jsonify({"error": "Topic not found"}), 404

@app.route('/api/class/<class_id>')
def get_class_data(class_id):
    return jsonify(data.get(class_id, {}))

@app.route('/api/classes')
def get_all_classes():
    return jsonify(list(data.keys()))

@app.route('/api/class/<class_id>/subjects')
def get_subjects(class_id):
    return jsonify(list(data.get(class_id, {}).keys()))

@app.route('/api/class/<class_id>/subject/<subject>/topics')
def get_topics(class_id, subject):
    subject_data = data.get(class_id, {}).get(subject, {})
    return jsonify(list(subject_data.keys()))

@app.route('/api/search/<query>')
def search_topics(query):
    results = []
    query_lower = query.lower()
    for class_id, subjects in data.items():
        for subject, topics in subjects.items():
            for topic, details in topics.items():
                if query_lower in topic.lower() or query_lower in details.get("description", "").lower():
                    results.append({
                        "class": class_id,
                        "subject": subject,
                        "topic": topic,
                        "difficulty": details.get("difficulty", ""),
                        "importance": details.get("importance", "")
                    })
    return jsonify(results)

@app.route('/api/difficulty/<level>')
def get_by_difficulty(level):
    results = []
    for class_id, subjects in data.items():
        for subject, topics in subjects.items():
            for topic, details in topics.items():
                if details.get("difficulty", "").lower() == level.lower():
                    results.append({
                        "class": class_id,
                        "subject": subject,
                        "topic": topic,
                        "importance": details.get("importance", "")
                    })
    return jsonify(results)

@app.route('/api/class/<class_id>/subject/<subject>/topic/<topic>/formulas')
def get_formulas(class_id, subject, topic):
    subject_data = data.get(class_id, {}).get(subject, {})
    for key in subject_data:
        if key.lower().replace(" ", "") == topic.lower().replace(" ", ""):
            return jsonify({"topic": key, "formulas": subject_data[key].get("formulas", [])})
    return jsonify({"error": "Topic not found"}), 404

@app.route('/api/class/<class_id>/subject/<subject>/topic/<topic>/roadmap')
def get_roadmap(class_id, subject, topic):
    subject_data = data.get(class_id, {}).get(subject, {})
    for key in subject_data:
        if key.lower().replace(" ", "") == topic.lower().replace(" ", ""):
            return jsonify({"topic": key, "roadmap": subject_data[key].get("roadmap", [])})
    return jsonify({"error": "Topic not found"}), 404

@app.route('/api/stats')
def get_stats():
    stats = {}
    for class_id, subjects in data.items():
        stats[class_id] = {}
        for subject, topics in subjects.items():
            difficulty_counts = {"Easy": 0, "Medium": 0, "Hard": 0}
            for topic, details in topics.items():
                d = details.get("difficulty", "")
                if d in difficulty_counts:
                    difficulty_counts[d] += 1
            stats[class_id][subject] = {
                "total_topics": len(topics),
                "difficulty_breakdown": difficulty_counts
            }
    return jsonify(stats)

# =========================================
# SHISHIMANU AI ASSISTANT — Gemini-Powered
# =========================================

@app.route('/api/shishimanu/chat', methods=['POST'])
def shishimanu_chat():
    if not _gemini_available:
        return jsonify({'error': 'google-generativeai package not installed. Run: pip install google-generativeai'}), 500

    api_key = os.environ.get('GEMINI_API_KEY')
    if not api_key:
        return jsonify({'error': 'GEMINI_API_KEY not set. Please set it as an environment variable.'}), 500

    body = request.get_json(force=True, silent=True) or {}
    user_message = body.get('message', '').strip()
    history = body.get('history', [])  # list of {role, content} dicts

    if not user_message:
        return jsonify({'error': 'No message provided'}), 400

    system_prompt = (
        "You are Shishimanu 🦁 — a friendly, enthusiastic, and brilliant AI study companion "
        "built into EDU-ARCHITECT, an educational platform for Class 9, 10, 11, and 12 students "
        "studying Physics, Chemistry, and Mathematics. "
        "Your personality is warm, encouraging, and a little playful (like your cute lion mascot!). "
        "You explain complex concepts in simple language with examples, help with formulas, "
        "solve problems step-by-step, and motivate students to learn. "
        "Keep responses concise and student-friendly. Use emojis occasionally to stay engaging. "
        "Always encourage curiosity and a growth mindset."
    )

    try:
        genai.configure(api_key=api_key)
        model = genai.GenerativeModel(
            model_name='gemini-1.5-flash',
            system_instruction=system_prompt
        )

        # Build Gemini chat history (all turns except the last user message)
        gemini_history = []
        for h in history[-10:]:
            role = h.get('role', 'user')
            content = h.get('content', '')
            if role == 'assistant':
                role = 'model'  # Gemini uses 'model' instead of 'assistant'
            if role in ('user', 'model') and content:
                gemini_history.append({'role': role, 'parts': [content]})

        chat = model.start_chat(history=gemini_history)
        response = chat.send_message(user_message)
        reply = response.text
        return jsonify({'reply': reply})
    except Exception as e:
        error_msg = str(e)
        if 'API_KEY_INVALID' in error_msg or 'invalid api key' in error_msg.lower():
            return jsonify({'error': 'Invalid API key. Please check your GEMINI_API_KEY.'}), 401
        return jsonify({'error': error_msg}), 500


@app.route('/favicon.ico')
def favicon():
    return '', 204

@app.route('/<path:path>')
def serve_static(path):
    return send_from_directory('.', path)

if __name__ == '__main__':
    app.run(debug=True, port=3000)
