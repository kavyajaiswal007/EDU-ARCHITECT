const CHAPTER_DATA = {
    "11": {
        "Physics": {
            "Physics & Measurement": {
                "difficulty": "Easy",
                "description": "Physics & Measurement is the gateway chapter of JEE Physics. It establishes the language of physics — units, dimensions, and the art of precise measurement. Every physical quantity can be expressed in terms of seven fundamental SI units. Dimensional analysis is a powerful tool to derive formulas, check equations, and convert units. Understanding significant figures and errors is essential for experimental accuracy, tested extensively in JEE practical-based questions.",
                "roadmap": ["Physical Quantities & SI Units", "Dimensional Analysis", "Significant Figures", "Errors in Measurement", "Instruments — Vernier Caliper & Screw Gauge"],
                "topics": ["Fundamental & Derived Quantities", "Dimensions of Physical Quantities", "Principle of Homogeneity", "Absolute, Relative & Percentage Error", "Propagation of Errors", "Least Count of Instruments", "Order of Magnitude"],
                "formulas": ["[M L T] notation for dimensions", "% Error = (ΔA/A) × 100", "If Z = A×B, ΔZ/Z = ΔA/A + ΔB/B", "Least count = MSD − VSD"],
                "applications": ["Checking correctness of equations", "Deriving unknown formulas", "Converting between unit systems"],
                "common_mistakes": ["Forgetting to square dimensions in area/volume", "Adding quantities with different dimensions", "Confusing precision with accuracy"],
                "importance": "5-7% JEE weightage. Easy marks — always appears in JEE Main. Master dimensional analysis."
            },
            "Vectors": {
                "difficulty": "Medium",
                "description": "Vectors are the mathematical backbone of all physics. Unlike scalars, vectors carry both magnitude and direction — making them essential for describing forces, velocities, and electric fields. This chapter covers vector addition, subtraction, dot product, cross product, and resolution. Mastery of vectors is non-negotiable for JEE since they appear in mechanics, electromagnetism, and optics.",
                "roadmap": ["Scalars vs Vectors", "Vector Addition (Triangle & Parallelogram Law)", "Resolution of Vectors", "Dot Product", "Cross Product"],
                "topics": ["Unit Vectors î, ĵ, k̂", "Vector Addition & Subtraction", "Component Method", "Scalar Product & properties", "Vector Product & Right-Hand Rule", "Position & Displacement Vectors"],
                "formulas": ["A·B = AB cosθ", "|A×B| = AB sinθ", "R² = A² + B² + 2AB cosθ", "tan α = B sinθ/(A + B cosθ)", "Magnitude = √(Ax² + Ay² + Az²)"],
                "applications": ["Finding resultant force", "Torque calculation (r × F)", "Work done (F · displacement)"],
                "common_mistakes": ["Confusing dot and cross product results", "Not using right-hand rule for cross product", "Adding magnitudes instead of vectors"],
                "importance": "Used in every chapter. 2-3 direct JEE questions. Essential throughout the syllabus."
            },
            "Motion in One Dimension": {
                "difficulty": "Easy",
                "description": "Motion in One Dimension is the study of objects moving along a straight line. It introduces displacement, velocity, and acceleration in a single dimension. The three equations of motion are the most-used formulas in all of physics. Velocity-time and displacement-time graphs are visual tools frequently tested in JEE Main.",
                "roadmap": ["Distance vs Displacement", "Speed vs Velocity", "Uniform & Non-uniform Acceleration", "Equations of Motion", "Graphs: s-t, v-t, a-t"],
                "topics": ["Average & Instantaneous Velocity", "Uniform Acceleration Equations", "Free Fall & Motion Under Gravity", "Relative Motion in 1D", "Analysis of v-t graphs", "Reaction Time"],
                "formulas": ["v = u + at", "s = ut + ½at²", "v² = u² + 2as", "s_nth = u + a(2n−1)/2", "g = 9.8 m/s²"],
                "applications": ["Free fall problems", "Particle motion from graphs", "Relative velocity between trains/cars"],
                "common_mistakes": ["Not taking sign convention consistently", "Confusing distance and displacement", "Misreading velocity-time graphs"],
                "importance": "Foundation of mechanics. 1-2 questions in JEE Main. Easy scoring with graph problems."
            },
            "Motion in Two Dimensions": {
                "difficulty": "Medium",
                "description": "Motion in Two Dimensions extends kinematics to a plane. The most important application is Projectile Motion — the path of a ball thrown at an angle. Uniform Circular Motion introduces centripetal acceleration and is the base for banking of roads and planetary orbits. Relative velocity in 2D (river-boat, rain-man problems) is a JEE favourite.",
                "roadmap": ["2D Kinematics — Component Method", "Projectile Motion", "Horizontal Projectile", "Uniform Circular Motion", "Relative Velocity in 2D"],
                "topics": ["Projectile Range, Height, Time of Flight", "Equation of Trajectory", "Centripetal & Centrifugal Acceleration", "Banking of Roads", "River-Boat Problems", "Conical Pendulum"],
                "formulas": ["R = u²sin2θ/g", "H = u²sin²θ/2g", "T = 2u sinθ/g", "aₖ = v²/r = ω²r", "tan θ = v²/rg (banking)"],
                "applications": ["Ball throw, cannon fire", "Satellite circular motion", "Car on banked curves"],
                "common_mistakes": ["Not resolving initial velocity into components", "Using g = 10 inconsistently", "Forgetting centripetal is always towards center"],
                "importance": "High weightage — 2-3 JEE questions. Projectile motion is one of the most tested topics."
            },
            "Laws of Motion": {
                "difficulty": "Hard",
                "description": "Newton's Laws of Motion are the cornerstone of classical mechanics. The First Law defines inertia. The Second Law (F = ma) quantifies force. The Third Law states action = reaction. This chapter covers pulley problems, connected bodies, friction on inclined planes, and pseudo forces in non-inertial frames — all JEE favourites.",
                "roadmap": ["Newton's Three Laws", "Free Body Diagram (FBD) Technique", "Friction — Static, Kinetic, Rolling", "Connected Bodies & Pulleys", "Pseudo Force in Non-Inertial Frames"],
                "topics": ["Inertia & Mass", "Normal Force", "Tension in Strings", "Static & Kinetic Friction", "Angle of Friction", "Atwood Machine", "Wedge Problems"],
                "formulas": ["F = ma", "p = mv", "Impulse J = FΔt = Δp", "f_s ≤ μ_s N", "f_k = μ_k N", "tan φ = μ"],
                "applications": ["Car braking distance", "Lift problems (apparent weight)", "Rocket propulsion"],
                "common_mistakes": ["Incomplete free body diagrams", "Applying pseudo force in inertial frame", "Wrong friction direction"],
                "importance": "15-20% of JEE mechanics questions. Mastering FBD is the single most important skill."
            },
            "Work, Energy & Power": {
                "difficulty": "Medium",
                "description": "Work, Energy & Power is one of the most elegant chapters. Work is done when a force produces displacement. The Work-Energy Theorem connects net work to change in KE. Conservation of mechanical energy bypasses force analysis entirely. Elastic and inelastic collisions and power round off this JEE-critical chapter.",
                "roadmap": ["Work by Constant & Variable Force", "Kinetic & Potential Energy", "Work-Energy Theorem", "Conservation of Mechanical Energy", "Power & Efficiency"],
                "topics": ["W = Fs cosθ", "Work done by spring force", "Conservative & Non-Conservative Forces", "Potential Energy Curve", "Elastic & Inelastic Collisions", "Coefficient of Restitution", "Instantaneous Power"],
                "formulas": ["W = F·s·cosθ", "W = ½kx² (spring)", "KE = ½mv²", "PE = mgh", "W_net = ΔKE", "P = F·v"],
                "applications": ["Roller coaster energy conservation", "Spring systems", "Car engine power output"],
                "common_mistakes": ["Using W = Fs without cosθ", "Not accounting for friction in energy conservation", "Confusing elastic vs inelastic collision"],
                "importance": "10-12% JEE weightage. Conservation of energy is the most powerful problem-solving tool."
            },
            "Centre of Mass, Momentum & Collisions": {
                "difficulty": "Hard",
                "description": "The Centre of Mass (COM) represents the entire mass distribution of a system. Conservation of linear momentum holds when no external force acts. Elastic collisions conserve both momentum and KE; inelastic do not. This chapter is essential for problems involving explosions, rocket propulsion, and oblique collisions.",
                "roadmap": ["Centre of Mass — Discrete & Continuous Systems", "Motion of COM", "Conservation of Linear Momentum", "Types of Collisions", "Coefficient of Restitution"],
                "topics": ["COM of triangles, semicircles, hemispheres", "Velocity & Acceleration of COM", "Explosion problems", "Perfectly Elastic Collision formulas", "Perfectly Inelastic Collision", "Rocket Propulsion"],
                "formulas": ["x_COM = Σmᵢxᵢ / Σmᵢ", "p_total = constant (no ext. force)", "e = 1 (elastic), e = 0 (perfectly inelastic)", "v_rocket = u ln(M₀/M)"],
                "applications": ["Explosion problems", "Ballistic pendulum", "Rocket thrust"],
                "common_mistakes": ["Applying momentum conservation when external forces exist", "Forgetting COM doesn't change during internal forces"],
                "importance": "8-10% JEE weightage. Collisions appear in almost every JEE paper."
            },
            "Rotational Motion": {
                "difficulty": "Hard",
                "description": "Rotational Motion is the most mathematically rich chapter in Class 11 Physics. Every translational concept has a rotational analogue: force → torque, mass → moment of inertia, linear momentum → angular momentum. The parallel and perpendicular axis theorems allow calculation of MI for complex bodies. Rolling motion combines translational and rotational KE.",
                "roadmap": ["Angular Kinematics (α, ω, θ)", "Torque & Couple", "Moment of Inertia & Theorems", "Angular Momentum & Conservation", "Rolling Motion Without Slipping"],
                "topics": ["Equations of Rotational Motion", "Torque = r × F", "MI of rod, disk, ring, sphere", "Parallel Axis Theorem: I = I_cm + Md²", "Perpendicular Axis Theorem: I_z = I_x + I_y", "Rolling condition v = Rω"],
                "formulas": ["τ = Iα", "L = Iω", "KE_roll = ½Iω² + ½mv²", "I_disk = ½MR²", "I_sphere = ⅖MR²", "I_ring = MR²", "I_rod = ML²/12"],
                "applications": ["Spinning top & gyroscope", "Figure skater pulling arms in", "Rolling objects down inclines"],
                "common_mistakes": ["Using Iα without rolling constraint", "Forgetting ½mv² for rolling KE", "Wrong axis for parallel axis theorem"],
                "importance": "15% JEE weightage. One of the highest-scoring chapters in JEE Advanced."
            },
            "Gravitation": {
                "difficulty": "Medium",
                "description": "Gravitation governs the motion of planets, moons, and satellites. Newton's Law gives the attractive force between any two masses. This chapter covers gravitational field, potential, escape velocity, orbital velocity, geostationary satellites, and Kepler's three laws. Variations of 'g' with altitude, depth, and latitude are tested frequently in JEE Main.",
                "roadmap": ["Newton's Law of Gravitation", "Gravitational Field & Potential", "Variation of g", "Orbital & Escape Velocity", "Satellites & Kepler's Laws"],
                "topics": ["g = F/m (field)", "V = −GM/r (potential)", "g' = g(1−2h/R) at height h", "g' = g(1−d/R) at depth d", "Orbital velocity v₀ = √(GM/r)", "Escape velocity vₑ = √(2GM/R)", "Geostationary Satellites (T = 24 hr)", "Kepler's 3 Laws"],
                "formulas": ["F = Gm₁m₂/r²", "g = GM/R²", "vₑ = √(2gR) = 11.2 km/s", "v₀ = √(gR)", "T² ∝ r³"],
                "applications": ["GPS & communication satellites", "Moon's orbital period", "Space mission trajectories"],
                "common_mistakes": ["Confusing g at surface vs at height h", "Using vₑ = √(2gR) only at Earth's surface"],
                "importance": "8-10% JEE weightage. Satellites and escape velocity appear every year."
            },
            "Mechanical Properties of Solids": {
                "difficulty": "Medium",
                "description": "Mechanical Properties of Solids explores how solid materials respond to external forces. Stress is the internal restoring force per unit area; strain is the fractional deformation. The three elastic moduli (Young's, Bulk, Shear) quantify material response. Hooke's Law states stress is proportional to strain within the elastic limit.",
                "roadmap": ["Stress & Strain", "Hooke's Law & Elastic Limit", "Young's Modulus", "Bulk & Shear Modulus", "Elastic Potential Energy"],
                "topics": ["Tensile, Compressive & Shear Stress", "Young's Modulus Y = Stress/Strain", "Bulk Modulus B = −P/(ΔV/V)", "Poisson's Ratio", "Elastic PE = ½ × Stress × Strain × Volume"],
                "formulas": ["Y = FL/AΔL", "B = −VΔP/ΔV", "Elastic PE = ½YAε²l"],
                "applications": ["Steel vs rubber elasticity", "Bridges and building materials", "Bone mechanics"],
                "common_mistakes": ["Confusing modulus of rigidity with Young's modulus", "Wrong elastic PE formula"],
                "importance": "3-5% JEE weightage. Short but conceptual. Know all 3 moduli."
            },
            "Mechanical Properties of Fluids": {
                "difficulty": "Medium",
                "description": "Mechanical Properties of Fluids covers pressure, buoyancy, and fluid dynamics. Archimedes' Principle explains buoyancy. Bernoulli's Theorem (conservation of energy in fluid flow) explains aircraft lift and venturimeter. Viscosity and surface tension govern everyday experiences from blood flow to soap bubbles.",
                "roadmap": ["Fluid Pressure & Pascal's Law", "Archimedes' Principle & Buoyancy", "Bernoulli's Theorem", "Viscosity & Stoke's Law", "Surface Tension & Capillarity"],
                "topics": ["P = ρgh", "Pascal's Law", "Buoyant Force = ρ_fluid × V_submerged × g", "Continuity: A₁v₁ = A₂v₂", "Torricelli's theorem", "Stoke's Law F = 6πηrv", "Terminal Velocity", "Excess Pressure in bubble: 4T/r", "Capillary Rise h = 2T cosθ/ρrg"],
                "formulas": ["P + ½ρv² + ρgh = const (Bernoulli)", "h = 2T cosθ/(ρrg)", "v_terminal = 2r²(ρ-σ)g/9η"],
                "applications": ["Aircraft wing lift", "Hydraulic brakes", "Blood flow through arteries"],
                "common_mistakes": ["Applying Bernoulli where flow is not streamline", "Confusing 2T vs 4T for drops vs bubbles"],
                "importance": "8-10% JEE weightage. Bernoulli and surface tension are frequently tested."
            },
            "Thermal Properties of Matter & Heat Transfer": {
                "difficulty": "Medium",
                "description": "This chapter covers heat, temperature scales, thermal expansion, specific heat, latent heat, and all three modes of heat transfer — conduction, convection, and radiation. Stefan-Boltzmann law, Wien's displacement law, and Newton's law of cooling are all tested in JEE.",
                "roadmap": ["Temperature Scales & Conversion", "Thermal Expansion", "Specific Heat & Calorimetry", "Latent Heat & Phase Changes", "Modes of Heat Transfer"],
                "topics": ["ΔL = αLΔT, ΔA = 2αAΔT, ΔV = 3αVΔT", "Q = mcΔT", "Q = mL", "Heat gained = Heat lost (calorimetry)", "Newton's Law of Cooling", "Stefan-Boltzmann E = σT⁴", "Wien's λ_max T = b"],
                "formulas": ["Q = mcΔT", "Q = mL", "dQ/dt = kA(ΔT/x)", "E = εσT⁴", "λ_max = b/T"],
                "applications": ["Thermos flask", "Land & sea breezes", "Pressure cooker"],
                "common_mistakes": ["Confusing specific heat with heat capacity", "Not converting to Kelvin for radiation"],
                "importance": "8-10% JEE weightage. Wien's and Stefan's laws appear regularly."
            },
            "Thermodynamics": {
                "difficulty": "Hard",
                "description": "Thermodynamics is the science of energy transformation. The Zeroth Law defines temperature; the First Law is conservation of energy; the Second Law introduces entropy and irreversibility. The Carnot engine gives maximum efficiency. PV diagrams for isothermal, adiabatic, isochoric, and isobaric processes are core JEE exam tools.",
                "roadmap": ["Thermodynamic Systems & Processes", "Zeroth & First Laws", "PV Diagrams — 4 Processes", "Second Law & Carnot Engine", "Entropy"],
                "topics": ["Internal Energy ΔU", "W = ∫PdV", "Isothermal: W = nRT ln(V₂/V₁)", "Adiabatic: PVγ = constant", "Isochoric: W = 0", "Isobaric: W = PΔV", "Carnot Cycle — 4 steps", "η = 1 − T_cold/T_hot"],
                "formulas": ["ΔU = Q − W", "W = PΔV (isobaric)", "W = nRT ln(V₂/V₁) (isothermal)", "PVγ = const (adiabatic)", "η_Carnot = 1 − T₂/T₁", "COP = T₂/(T₁−T₂)"],
                "applications": ["Car engines", "Refrigerators & ACs", "Steam turbines"],
                "common_mistakes": ["Wrong sign convention for Q and W", "Using adiabatic formula for isothermal"],
                "importance": "10-12% JEE weightage. Carnot efficiency and PV diagrams appear in almost every JEE paper."
            },
            "Kinetic Theory of Gases": {
                "difficulty": "Medium",
                "description": "The Kinetic Theory provides the microscopic explanation for macroscopic gas behaviour. It models gas as tiny particles in random elastic collisions. From this model we derive pressure, temperature, the ideal gas equation PV = nRT, and Maxwell speed distributions. Degrees of freedom and equipartition theorem explain heat capacities.",
                "roadmap": ["Postulates of Kinetic Theory", "Pressure & Temperature from Molecular Model", "Gas Laws", "RMS, Average & Most Probable Speed", "Degrees of Freedom & Equipartition"],
                "topics": ["PV = nRT", "KE = 3/2 kT per molecule", "v_rms = √(3RT/M)", "v_avg = √(8RT/πM)", "v_mp = √(2RT/M)", "Mean Free Path", "Cv = fR/2, Cp = (f+2)R/2, γ = Cp/Cv"],
                "formulas": ["PV = nRT", "KE = 3/2 kT", "v_rms : v_avg : v_mp = √3 : √(8/π) : √2", "γ = 5/3 (monoatomic), 7/5 (diatomic)"],
                "applications": ["Explaining gas laws from molecular model", "Gas thermometers", "Effusion of gases"],
                "common_mistakes": ["Confusing v_rms, v_avg, and v_mp", "Wrong degrees of freedom for diatomic vs monoatomic"],
                "importance": "8% JEE weightage. Speed ratios and Cv, Cp, γ are must-know."
            },
            "Waves & Sound": {
                "difficulty": "Medium",
                "description": "Waves chapter explores how energy propagates through a medium. Waves are transverse (light) or longitudinal (sound). The wave equation y = A sin(kx − ωt) encodes amplitude, frequency, and speed. Superposition creates interference, beats, and standing waves. Doppler effect explains the change in frequency with relative motion.",
                "roadmap": ["Wave Parameters", "Wave Equation & Phase", "Superposition & Interference", "Standing Waves — Strings & Pipes", "Beats & Doppler Effect"],
                "topics": ["v = fλ", "Phase difference: Δφ = (2π/λ)Δx", "Nodes & Antinodes", "Harmonics & Overtones", "Open & Closed Pipe", "Beat frequency = |f₁ − f₂|", "Doppler: f' = f(v±v_observer)/(v∓v_source)"],
                "formulas": ["v = fλ", "y = A sin(kx − ωt)", "v_sound = √(γP/ρ)", "f_n (string) = n/2L √(T/μ)", "f_n (open pipe) = nv/2L", "f_n (closed pipe) = (2n−1)v/4L"],
                "applications": ["Musical instruments", "Medical ultrasound", "Speed gun (Doppler radar)"],
                "common_mistakes": ["Forgetting closed pipe has only odd harmonics", "Wrong sign in Doppler formula"],
                "importance": "8-10% JEE weightage. Standing waves and Doppler are annual JEE favourites."
            },
            "Simple Harmonic Motion": {
                "difficulty": "Medium",
                "description": "Simple Harmonic Motion (SHM) occurs when a restoring force acts proportional to displacement: F = −kx. Motion is sinusoidal with amplitude, time period, frequency, and phase. Spring-mass systems and simple pendulums are key oscillators. Energy in SHM continuously exchanges between KE and PE, with total energy conserved.",
                "roadmap": ["Conditions for SHM", "Displacement, Velocity, Acceleration in SHM", "Spring-Mass System", "Simple Pendulum", "Energy in SHM"],
                "topics": ["x = A sin(ωt + φ)", "v = Aω cos(ωt + φ)", "a = −ω²x", "T = 2π√(m/k)", "T = 2π√(l/g)", "KE = ½mω²(A²−x²)", "PE = ½mω²x²", "Damped & Forced Oscillations"],
                "formulas": ["x = A sin(ωt + φ)", "T = 2π/ω", "T_spring = 2π√(m/k)", "T_pendulum = 2π√(L/g)", "E_total = ½kA²"],
                "applications": ["Pendulum clock", "Car suspension", "Resonance bridges"],
                "common_mistakes": ["Confusing time period with frequency", "Forgetting amplitude is max displacement from equilibrium"],
                "importance": "10-12% JEE weightage. SHM combined with waves and LC circuits in advanced problems."
            }
        },

        "Chemistry": {
            "Some Basic Concepts of Chemistry": {
                "difficulty": "Easy",
                "description": "Some Basic Concepts of Chemistry lays the quantitative foundation of chemistry. It covers Dalton's Atomic Theory, laws of chemical combination, SI units, mole concept, molar mass, and stoichiometry. The mole bridges the atomic world (atoms, molecules) with the macroscopic world (grams, litres).",
                "roadmap": ["Matter & its Nature", "Laws of Chemical Combination", "Atomic & Molecular Masses", "Mole Concept & Molar Mass", "Stoichiometry & Limiting Reagent"],
                "topics": ["Dalton's Atomic Theory", "Law of Conservation of Mass", "Law of Definite Proportions", "Law of Multiple Proportions", "Avogadro's Law", "Mole = 6.022×10²³ particles", "Molar mass", "Percentage composition", "Empirical & Molecular formula", "Limiting reagent", "% yield"],
                "formulas": ["Moles = mass/molar mass", "Moles = volume(L)/22.4 (at STP)", "% composition = (mass of element / molar mass) × 100"],
                "applications": ["Titration calculations", "Finding limiting reagent", "Industrial chemical production"],
                "common_mistakes": ["Confusing empirical and molecular formula", "Forgetting to convert g to mol before stoichiometry", "Wrong significant figures"],
                "importance": "8-10% JEE Main weightage. Easy marks — mole concept and stoichiometry appear every year."
            },
            "Structure of Atom": {
                "difficulty": "Hard",
                "description": "Structure of Atom traces the development of atomic models from Thomson and Rutherford to Bohr and the quantum mechanical model. Quantum numbers (n, l, ml, ms) and orbital shapes (s, p, d) define each electron's state and determine chemical bonding.",
                "roadmap": ["Thomson & Rutherford Models", "Bohr's Model of Hydrogen", "Dual Nature & Uncertainty Principle", "Quantum Mechanical Model", "Electronic Configuration"],
                "topics": ["Rutherford's gold foil experiment", "Bohr's postulates: rn = 0.529 n² Å, En = -13.6/n² eV", "Hydrogen spectrum — Lyman, Balmer, Paschen series", "de Broglie: λ = h/mv", "Heisenberg: Δx·Δp ≥ h/4π", "Quantum numbers: n, l, ml, ms", "Shapes of s, p, d orbitals", "Aufbau principle, Pauli's exclusion, Hund's rule", "Extra stability of half-filled and fully filled orbitals"],
                "formulas": ["E = hν", "λ = h/p = h/mv (de Broglie)", "Δx·Δp ≥ h/4π (Heisenberg)", "1/λ = R(1/n₁² - 1/n₂²) (Rydberg)", "En = -13.6/n² eV (Bohr, H atom)"],
                "applications": ["Flame tests (atomic spectra)", "Laser technology", "Electron microscopes"],
                "common_mistakes": ["Forgetting max electrons per orbital = 2", "Wrong orbital filling order (remember: 4s before 3d)", "Confusing quantum numbers l and ml"],
                "importance": "10-12% JEE weightage. Electronic config and quantum numbers appear in every JEE paper."
            },
            "Classification of Elements & Periodicity": {
                "difficulty": "Medium",
                "description": "Classification of Elements & Periodicity explains the logic behind the periodic table — Mendeleev's arrangement and Moseley's atomic number discovery. Periodic trends (atomic radius, ionization energy, electron affinity, electronegativity) predict reactivity patterns across groups and periods.",
                "roadmap": ["Historical Development of Periodic Table", "Modern Periodic Law", "Blocks (s, p, d, f)", "Periodic Trends — Atomic Radius", "Periodic Trends — IE, EA, EN"],
                "topics": ["Modern periodic law: properties are periodic functions of atomic number", "s, p, d, f block elements", "Atomic radius trends: increases down group, decreases across period", "Ionic radius: cation < parent atom < anion", "Ionization energy: IE₁ < IE₂ < IE₃", "Electronegativity: Pauling scale", "Metallic & non-metallic character", "Diagonal relationship (Li-Mg, Be-Al)"],
                "formulas": ["Effective nuclear charge (Zeff) = Z - σ (Slater's rules)"],
                "applications": ["Predicting compound formation", "Explaining anomalous first ionization energies", "Diagonal relationships in industry"],
                "common_mistakes": ["Forgetting EA is exothermic (negative)", "Wrong trend for electron affinity (not always regular)", "Confusing atomic and ionic radius trends"],
                "importance": "8-10% JEE weightage. Mostly conceptual — easy marks with good memory."
            },
            "Chemical Bonding & Molecular Structure": {
                "difficulty": "Hard",
                "description": "Chemical Bonding & Molecular Structure explains why and how atoms bond. VSEPR theory predicts molecular geometry. Hybridization (sp, sp², sp³, sp³d, sp³d²) determines molecular shape. MO Theory explains paramagnetism and bond order. Hydrogen bonding, Fajan's rules, and dipole moment are essential JEE topics.",
                "roadmap": ["Octet Rule & Exceptions", "Ionic vs Covalent Bonding", "Lewis Dot Structures", "VSEPR Theory — Molecular Geometry", "Hybridization & MO Theory"],
                "topics": ["Electronegativity & Fajan's rules", "Lewis dot structures", "Formal charge: FC = V - L - B/2", "VSEPR — shapes: linear, bent, trigonal planar, tetrahedral, trigonal bipyramidal, octahedral", "Hybridization: sp (linear), sp² (planar), sp³ (tetrahedral)", "Bond order = (Nb - Na)/2", "Paramagnetism of O₂ (explained by MO)", "Hydrogen bonding — intra and inter molecular", "Dipole moment"],
                "formulas": ["Bond order = (Nb - Na)/2", "FC = V - L - B/2", "% Ionic character = 16|Δχ| + 3.5(Δχ)²"],
                "applications": ["Predicting molecular shapes", "Explaining boiling points via H-bonding", "DNA double helix (H-bonds)"],
                "common_mistakes": ["Wrong number of lone pairs in VSEPR", "Confusing sigma/pi with single/double bonds", "Forgetting O₂ is paramagnetic"],
                "importance": "12-15% JEE weightage. VSEPR and MO theory problems appear every year."
            },
            "States of Matter": {
                "difficulty": "Medium",
                "description": "States of Matter explores how matter behaves as a gas, liquid, or solid. Gas laws (Boyle's, Charles's, Graham's, Dalton's, Ideal Gas Law PV = nRT), kinetic molecular theory, and van der Waals equation explain gas behavior. Surface tension, viscosity, and vapour pressure describe the liquid state.",
                "roadmap": ["Gas Laws — Boyle, Charles, Graham", "Ideal Gas Equation", "Dalton's Law of Partial Pressure", "Kinetic Theory of Gases", "Real Gases & van der Waals"],
                "topics": ["Boyle's Law: PV = constant (constant T)", "Charles's Law: V/T = constant (constant P)", "Avogadro's Law: V ∝ n", "Graham's Law: r₁/r₂ = √(M₂/M₁)", "Dalton's Law: P_total = P₁ + P₂ + P₃...", "PV = nRT (Ideal Gas)", "Compressibility factor Z = PV/nRT", "van der Waals equation: (P + an²/V²)(V - nb) = nRT", "Surface tension, viscosity, vapour pressure"],
                "formulas": ["PV = nRT", "r₁/r₂ = √(M₂/M₁)", "(P + an²/V²)(V - nb) = nRT", "Z = PV/nRT", "vrms = √(3RT/M)"],
                "applications": ["Weather balloons", "Scuba diving (Dalton's Law)", "Industrial gas storage"],
                "common_mistakes": ["Using wrong R value (8.314 vs 0.0821 L·atm/mol·K)", "Forgetting to convert temperature to Kelvin"],
                "importance": "8-10% JEE weightage. Gas law numericals are straightforward — easy marks."
            },
            "Thermodynamics (Chemistry)": {
                "difficulty": "Hard",
                "description": "Chemical Thermodynamics studies energy changes in chemical processes. First Law introduces internal energy (U), heat (q), and work (w). Hess's Law allows indirect enthalpy calculations. Gibbs Free Energy (G = H - TS) determines spontaneity — a central concept for chemical equilibrium.",
                "roadmap": ["System, Surroundings & State Functions", "First Law — Internal Energy & Enthalpy", "Hess's Law", "Second Law — Entropy", "Gibbs Free Energy & Spontaneity"],
                "topics": ["Open, closed, isolated systems", "State functions: U, H, S, G", "First Law: ΔU = q + w", "Work: w = -PΔV (expansion)", "Enthalpy: ΔH = ΔU + PΔV", "Exothermic (ΔH < 0) & Endothermic (ΔH > 0)", "Hess's Law: ΔHrxn = ΣΔHproducts - ΣΔHreactants", "Gibbs: ΔG = ΔH - TΔS", "ΔG° = -RT ln K = -nFE°"],
                "formulas": ["ΔU = q + w", "ΔH = ΔU + ΔnRT (for gases)", "ΔG = ΔH - TΔS", "ΔG° = -RT ln K", "ΔG° = -nFE°"],
                "applications": ["Predicting reaction spontaneity", "Fuel cell efficiency", "Industrial process design"],
                "common_mistakes": ["Sign error: work done ON system is +w", "Forgetting ΔG < 0 for spontaneous reaction", "Wrong sign for enthalpy of combustion"],
                "importance": "12-15% JEE weightage. High scoring in JEE Advanced — ΔG = -RT ln K is critical."
            },
            "Equilibrium": {
                "difficulty": "Hard",
                "description": "Equilibrium is one of the most concept-heavy JEE Chemistry chapters. Kc and Kp quantify the equilibrium position. Le Chatelier's Principle predicts equilibrium shifts. Ionic equilibrium covers weak acid/base dissociation, pH, buffer solutions (Henderson-Hasselbalch), solubility product Ksp, and common ion effect.",
                "roadmap": ["Law of Chemical Equilibrium", "Kc & Kp", "Le Chatelier's Principle", "Ionic Equilibrium — Acids & Bases", "Buffer & Solubility"],
                "topics": ["Dynamic equilibrium", "Kp = Kc(RT)^Δn", "Reaction quotient Q: Q < K (forward), Q > K (reverse)", "Le Chatelier: effect of concentration, pressure, temperature, catalyst", "Brønsted-Lowry & Lewis acids/bases", "Conjugate acid-base pairs", "Ka, Kb, Kw = Ka × Kb = 10⁻¹⁴", "pH = -log[H⁺]", "Buffer: pH = pKa + log([A⁻]/[HA])", "Ksp and solubility", "Common ion effect"],
                "formulas": ["Kp = Kc(RT)^Δn", "pH = -log[H⁺]", "Kw = [H⁺][OH⁻] = 10⁻¹⁴", "Henderson-Hasselbalch: pH = pKa + log([A⁻]/[HA])", "Ksp = [A⁺][B⁻] for AB"],
                "applications": ["Industrial Haber process (NH₃ synthesis)", "Buffer in blood (carbonate buffer pH 7.4)", "Stalagmites and stalactites (Ksp)"],
                "common_mistakes": ["Including solids/liquids in Kc/Kp expressions", "Confusing Q with K", "Wrong buffer pH formula"],
                "importance": "15-18% JEE weightage. Maximum JEE questions come from this chapter."
            },
            "Redox Reactions": {
                "difficulty": "Medium",
                "description": "Redox Reactions involve electron transfer between chemical species. Oxidation is loss of electrons (OIL), reduction is gain (RIG). Oxidation numbers allow identification of redox changes and balancing complex equations using the half-reaction method. Essential foundation for Class 12 electrochemistry.",
                "roadmap": ["Oxidation & Reduction Definitions", "Oxidation Number Rules", "Identifying Redox Reactions", "Balancing Redox Equations — Half-Reaction Method", "Disproportionation Reactions"],
                "topics": ["OIL-RIG mnemonic", "Oxidation number rules (O = -2, H = +1 etc.)", "Half-reaction method in acidic/basic medium", "Disproportionation: same element oxidized & reduced", "Common oxidizing agents: KMnO₄, K₂Cr₂O₇, H₂O₂, HNO₃", "Common reducing agents: Fe, C, Na₂S₂O₃"],
                "formulas": ["Electrons lost = electrons gained", "Net ionic equation method"],
                "applications": ["Corrosion of metals (rusting)", "Bleaching powder action", "Metabolic reactions (respiration)"],
                "common_mistakes": ["Wrong oxidation number for S, N, Cl (variable valency)", "Forgetting to balance charge with H⁺ or OH⁻", "Misidentifying oxidizing vs reducing agent"],
                "importance": "8-10% JEE weightage. Essential foundation for electrochemistry in Class 12."
            },
            "Hydrogen": {
                "difficulty": "Easy",
                "description": "Hydrogen is the simplest and most abundant element. This chapter covers hydrogen's unique position in the periodic table (resembles both alkali metals and halogens), isotopes (protium, deuterium, tritium), properties of water and heavy water (D₂O), and hydrogen peroxide (H₂O₂) as a mild oxidizer.",
                "roadmap": ["Position of H in Periodic Table", "Isotopes of Hydrogen", "Preparation & Properties of H₂", "Water & Heavy Water", "Hydrogen Peroxide (H₂O₂)"],
                "topics": ["Electronic configuration 1s¹ — resembles IA and VIIA", "Protium (H), Deuterium (D), Tritium (T)", "Water: bent shape (V-shaped), bond angle 104.5°, H-bonding", "Hard & Soft water, water treatment", "H₂O₂: mild oxidizing and reducing agent, bleaching", "H₂O₂ structure: non-planar, O-O single bond", "Hydrogen as fuel — fuel cells"],
                "formulas": ["H₂O₂ volume strength = molarity × 11.2"],
                "applications": ["Hydrogen fuel cells for clean energy", "Heavy water in nuclear reactors", "H₂O₂ in antiseptics and bleaching"],
                "common_mistakes": ["Forgetting H₂O₂ is a covalent compound (not ionic)", "Wrong structure of H₂O₂ (non-planar, not planar)"],
                "importance": "5-7% JEE weightage. Short chapter — factual questions in JEE Main."
            },
            "s-Block Elements (Alkali & Alkaline Earth Metals)": {
                "difficulty": "Medium",
                "description": "s-Block elements cover Groups 1 (Li, Na, K, Rb, Cs) and Group 2 (Be, Mg, Ca, Sr, Ba). The anomalous behavior of Li (resembles Mg — diagonal relationship) and Be (resembles Al) are important. Compounds NaOH, Na₂CO₃, NaCl, CaO, CaCO₃, and plaster of paris are tested extensively.",
                "roadmap": ["Group 1 — Alkali Metals", "Group 2 — Alkaline Earth Metals", "Anomalous Behavior of Li & Be", "Important Compounds", "Biological Significance"],
                "topics": ["ns¹ (Group 1) and ns² (Group 2) configurations", "Li anomalous: small size, resembles Mg (diagonal relationship)", "Be anomalous: resembles Al (diagonal relationship)", "Flame colors: Li-red, Na-yellow, K-violet, Ca-brick red, Ba-apple green", "NaOH (caustic soda) — chlor-alkali process", "Na₂CO₃ (washing soda), NaHCO₃ (baking soda)", "CaO (quicklime), Ca(OH)₂ (slaked lime), CaCO₃ (limestone)", "Plaster of Paris: 2CaSO₄·H₂O — sets with water"],
                "formulas": ["2Na + 2H₂O → 2NaOH + H₂", "CaO + H₂O → Ca(OH)₂ (slaking)"],
                "applications": ["Na as heat exchanger in nuclear reactors", "Mg in aircraft alloys", "CaCO₃ in cement production"],
                "common_mistakes": ["Confusing Na₂CO₃ (washing soda) with NaHCO₃ (baking soda)", "Forgetting Be is amphoteric (not basic like other Group 2)"],
                "importance": "8-10% JEE weightage. Important compound properties are frequently tested."
            },
            "p-Block Elements (Groups 13-18)": {
                "difficulty": "Hard",
                "description": "p-Block Elements spans Groups 13 to 18 — from reactive non-metals to inert noble gases. Covers boron (diborane), aluminium, carbon allotropes (Group 14), nitrogen compounds: NH₃ (Haber), HNO₃ (Ostwald), sulphuric acid (Contact process), halogens, and noble gas compounds (XeF₂, XeF₄).",
                "roadmap": ["Group 13 — B & Al compounds", "Group 14 — C & Si", "Group 15 — N, P, NH₃, HNO₃", "Group 16 — O, S, H₂SO₄", "Group 17 — Halogens & Group 18 — Noble Gases"],
                "topics": ["Boron: diborane (B₂H₆), boric acid (H₃BO₃)", "Carbon allotropes: diamond, graphite, fullerene", "NH₃: preparation (Haber process), pyramidal structure, basic", "HNO₃: preparation (Ostwald process), both oxidizer and acid", "PCl₃, PCl₅: structures and hydrolysis", "SO₂, SO₃, H₂SO₄ — contact process", "Ozone (O₃): bent, oxidizing, UV shield", "Halogens: F₂ strongest oxidizer, HF — weak acid", "Interhalogen compounds: ClF₃, BrF₅, IF₇", "Noble gases: Xe compounds — XeF₂, XeF₄, XeF₆"],
                "formulas": ["Haber: N₂ + 3H₂ → 2NH₃ (200 atm, 500°C, Fe catalyst)", "Ostwald: 4NH₃ + 5O₂ → 4NO + 6H₂O (Pt catalyst)", "Contact: 2SO₂ + O₂ → 2SO₃ (V₂O₅ catalyst)"],
                "applications": ["H₂SO₄ — 'king of chemicals' in fertilizer industry", "NH₃ — fertilizers, explosives", "Cl₂ — water purification"],
                "common_mistakes": ["Wrong hybridization for PCl₅ (sp³d) vs PH₃ (sp³)", "Forgetting HF is a weak acid despite F being most electronegative"],
                "importance": "15-18% JEE weightage. Largest chapter — most questions in JEE Main inorganic."
            },
            "Organic Chemistry — Basic Principles & Techniques": {
                "difficulty": "Medium",
                "description": "This chapter is the entry point into organic chemistry. Covers IUPAC nomenclature, types of isomerism (structural and stereoisomerism), reaction intermediates (carbocations, carbanions, free radicals), electronic effects (inductive, mesomeric, hyperconjugation), and basic types of organic reactions.",
                "roadmap": ["IUPAC Nomenclature", "Structural & Stereoisomerism", "Reaction Intermediates", "Electronic Effects", "Types of Organic Reactions"],
                "topics": ["IUPAC naming of alkanes, alkenes, alkynes, alcohols, aldehydes, acids", "Chain, position, functional group, metamerism isomers", "Geometric (cis-trans) isomerism: different groups on each C", "Optical isomerism: chiral centers, R & S configuration", "Carbocations: stability 3° > 2° > 1°", "Carbanions: stability 1° > 2° > 3°", "Free radicals: stability 3° > 2° > 1°", "Inductive effect: +I and -I groups", "Resonance/mesomeric effect: +M and -M", "Hyperconjugation — stability of alkenes and carbocations"],
                "formulas": ["+I effect: electron releasing groups", "-I effect: electron withdrawing groups", "+M: lone pairs donate into π system", "-M: withdraw electrons from π system"],
                "applications": ["Drug design (stereoisomers — thalidomide tragedy)", "Reaction mechanism prediction", "Industrial polymer design"],
                "common_mistakes": ["Wrong IUPAC name — not finding longest chain with highest substituents", "Confusing inductive with mesomeric effect", "Wrong stability order of carbocations"],
                "importance": "10-12% JEE weightage. Foundation for entire organic chemistry — master this first."
            },
            "Hydrocarbons": {
                "difficulty": "Hard",
                "description": "Hydrocarbons are organic compounds containing only C and H: alkanes (CnH₂n₊₂), alkenes (CnH₂n), alkynes (CnH₂n₋₂), and arenes (benzene). Each class has characteristic reactions: free radical substitution (alkanes), electrophilic addition following Markovnikov's rule (alkenes), and electrophilic aromatic substitution (benzene).",
                "roadmap": ["Alkanes — Preparation & Free Radical Substitution", "Alkenes — Electrophilic Addition & Markovnikov", "Alkynes — Acidic Character & Addition", "Benzene — Aromaticity & EAS", "Conformations & Cycloalkanes"],
                "topics": ["Alkanes: CnH₂n₊₂, sp³ C, chair & boat cyclohexane", "Preparation: Wurtz reaction, Kolbe's electrolysis", "Free radical halogenation: Cl₂/hv, reactivity 3° > 2° > 1°", "Markovnikov's Rule: H adds to more H-bearing C", "Ozonolysis, hydrogenation, hydration, polymerization", "Alkynes: sp C, acidic H in terminal alkynes, form metal acetylides", "Lindlar's catalyst: partial hydrogenation → cis alkene", "Benzene: 3 double bonds, resonance, Hückel's rule (4n+2 π electrons)", "EAS: nitration, halogenation, sulfonation, Friedel-Crafts reactions", "Ortho-para directors (+M groups) vs meta directors (-M groups)"],
                "formulas": ["Hückel: 4n+2 π electrons for aromaticity", "Markovnikov's rule", "Ozonolysis: R₁CH=CHR₂ + O₃ → R₁CHO + R₂CHO"],
                "applications": ["Petrol/diesel (alkanes)", "Ethylene oxide in antifreeze", "Benzene in pharmaceuticals and dyes"],
                "common_mistakes": ["Anti-Markovnikov product with HBr + peroxide (free radical)", "Wrong product in ozonolysis", "Forgetting benzene's EAS (not addition) reactions"],
                "importance": "15-18% JEE weightage. Hydrocarbons reactions are tested heavily every year."
            },
            "Environmental Chemistry": {
                "difficulty": "Easy",
                "description": "Environmental Chemistry deals with chemical processes in the environment and effects of human activities on air, water, and soil quality. Covers tropospheric/stratospheric pollution, greenhouse effect, ozone layer depletion by CFCs, acid rain, water pollutants (BOD), and eutrophication.",
                "roadmap": ["Atmospheric Pollution", "Greenhouse Effect & Global Warming", "Ozone Layer & Depletion", "Water Pollution", "Soil Pollution & Strategies"],
                "topics": ["Pollutants: primary (CO, SO₂, NO₂) and secondary (O₃, PAN)", "CO: toxic, binds haemoglobin, forms HbCO", "Greenhouse gases: CO₂, CH₄, N₂O, CFCs — trap IR radiation", "Ozone: O₃ in stratosphere blocks UV-B", "Ozone depletion by CFCs: Cl· radical chain reaction", "Acid rain: H₂SO₄ + HNO₃ from SO₂ + NO₂", "BOD (Biochemical Oxygen Demand) — water quality measure", "Eutrophication: algal bloom from excess N and P"],
                "formulas": ["BOD = dissolved O₂ consumed by bacteria in 5 days at 20°C"],
                "applications": ["CFC phase-out (Montreal Protocol)", "Catalytic converters in cars", "Sewage treatment plants"],
                "common_mistakes": ["Confusing greenhouse effect (thermal blanket) with ozone depletion (UV)", "Not knowing BOD is inversely related to water quality"],
                "importance": "3-5% JEE weightage. Short conceptual chapter — easy marks in JEE Main."
            }
        },

        "Maths": {
            "Trigonometry": {
                "difficulty": "Medium",
                "description": "Trigonometry covers angle measurement, trigonometric ratios, identities and equations, and inverse functions essential for calculus and problem-solving.",
                "roadmap": ["Angle Measurement", "Trigonometric Ratios", "Identities & Equations", "Inverse Functions", "Heights & Distances"],
                "topics": ["Radian & Degree Conversion", "Compound Angles", "Multiple & Sub-multiple Angles", "Trigonometric Equations", "Inverse Trig Functions"],
                "formulas": ["sin²θ + cos²θ = 1", "1 + tan²θ = sec²θ", "sin2θ = 2sinθcosθ", "cos2θ = cos²θ - sin²θ"],
                "importance": "Base for calculus. 12-14% board weightage. Essential JEE topic."
            },
            "Sets & Relations": {
                "difficulty": "Easy",
                "description": "Sets & Relations covers set theory basics, operations on sets, types of relations and functions, essential for building mathematical foundations.",
                "roadmap": ["Set Theory Basics", "Operations on Sets", "Relations", "Types of Relations", "Functions"],
                "topics": ["Types of Sets", "Venn Diagrams", "Cartesian Product", "Domain & Range", "Types of Functions"],
                "formulas": ["n(A∪B) = n(A) + n(B) - n(A∩B)", "n(A×B) = n(A) × n(B)"],
                "importance": "Logical foundation. 8% board weightage. Easy scoring."
            },
            "Complex Numbers": {
                "difficulty": "Hard",
                "description": "Complex Numbers extends the real number system. Covers imaginary numbers, algebra of complex numbers, Argand plane, polar form, and De Moivre's theorem.",
                "roadmap": ["Imaginary Numbers", "Algebra of Complex Numbers", "Argand Plane", "Polar Form", "De Moivre's Theorem"],
                "topics": ["i = √(-1)", "Operations on Complex Numbers", "Modulus & Argument", "Euler's Formula", "Roots of Unity"],
                "formulas": ["i² = -1", "|z| = √(a² + b²)", "z = r(cosθ + isinθ)", "(cosθ + isinθ)ⁿ = cosnθ + isinnθ"],
                "importance": "10-12% board weightage. High-value JEE problems."
            },
            "Permutations & Combinations": {
                "difficulty": "Medium",
                "description": "Permutations & Combinations covers counting techniques — the fundamental principle, permutations, combinations, and circular arrangements.",
                "roadmap": ["Fundamental Principle", "Permutations", "Combinations", "Circular Arrangements", "Restricted Cases"],
                "topics": ["Multiplication & Addition Principle", "nPr", "nCr", "Circular Permutations", "Selection with Constraints"],
                "formulas": ["nPr = n!/(n-r)!", "nCr = n!/[r!(n-r)!]", "nCr = nCn-r"],
                "importance": "10% board weightage. Probability base. Appears in JEE."
            },
            "Binomial Theorem": {
                "difficulty": "Medium",
                "description": "The Binomial Theorem gives a formula to expand (a+b)ⁿ. Covers expansion, general and middle terms, binomial coefficients, and Pascal's Triangle.",
                "roadmap": ["Expansion Formula", "General & Middle Terms", "Binomial Coefficients", "Applications"],
                "topics": ["(a+b)ⁿ Expansion", "Tr+1 Term", "Greatest Coefficient", "Pascal's Triangle"],
                "formulas": ["(a+b)ⁿ = Σ nCr·aⁿ⁻ʳ·bʳ", "Tr+1 = nCr·aⁿ⁻ʳ·bʳ"],
                "importance": "8-10% board weightage. Simple mark-scoring chapter."
            },
            "Sequences & Series": {
                "difficulty": "Hard",
                "description": "Sequences & Series covers AP, GP, HP, and special series. Includes nth term formulas, sum of n terms, infinite GP, and AM-GM-HM inequalities.",
                "roadmap": ["AP - Arithmetic Progression", "GP - Geometric Progression", "HP - Harmonic Progression", "Special Series", "AM, GM, HM"],
                "topics": ["nth Term Formulas", "Sum of n Terms", "Infinite GP", "Arithmetic Mean", "Geometric Mean Inequality"],
                "formulas": ["Tn = a + (n-1)d", "Sn = n/2[2a + (n-1)d]", "Tn = arⁿ⁻¹", "S∞ = a/(1-r) for |r|<1", "GM² = AM × HM"],
                "importance": "12% board weightage. Crucial for series convergence in calculus."
            },
            "Straight Lines": {
                "difficulty": "Medium",
                "description": "Straight Lines covers coordinate geometry — slope, various forms of line equations, angle between lines, and distance from a point to a line.",
                "roadmap": ["Distance & Section Formula", "Slope of Line", "Various Forms of Line", "Angle Between Lines", "Distance from Point to Line"],
                "topics": ["Slope-Intercept Form", "Point-Slope Form", "Two-Point Form", "Parallel & Perpendicular Lines", "Foot of Perpendicular"],
                "formulas": ["m = (y₂-y₁)/(x₂-x₁)", "y = mx + c", "ax + by + c = 0", "d = |ax₁+by₁+c|/√(a²+b²)"],
                "importance": "10% board weightage. Geometry foundation for JEE."
            },
            "Conic Sections": {
                "difficulty": "Hard",
                "description": "Conic Sections covers circle, parabola, ellipse, and hyperbola — standard forms, focus, directrix, eccentricity, tangent and normal.",
                "roadmap": ["Circle", "Parabola", "Ellipse", "Hyperbola", "General Equation"],
                "topics": ["Standard Forms", "Focus & Directrix", "Eccentricity", "Latus Rectum", "Tangent & Normal"],
                "formulas": ["x² + y² = r²", "y² = 4ax", "x²/a² + y²/b² = 1", "x²/a² - y²/b² = 1"],
                "importance": "15% board weightage. Maximum JEE questions from coordinate geometry."
            },
            "Limits & Derivatives": {
                "difficulty": "Hard",
                "description": "Limits & Derivatives covers the concept of limits, evaluation techniques, continuity, the derivative definition, and differentiation rules.",
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
                "description": "Electric Charges & Fields is the first chapter of Class 12 Physics and the foundation of electrostatics. It introduces the concept of electric charge, Coulomb's law for force between charges, and the electric field — a vector field that represents the force per unit charge at any point in space. Gauss's Law relates the electric flux through a closed surface to the enclosed charge, simplifying field calculations for symmetric configurations.",
                "roadmap": ["Properties of Electric Charge", "Coulomb's Law", "Electric Field & Field Lines", "Electric Dipole", "Gauss's Law & Applications"],
                "topics": ["Charge Quantization: q = ne", "Superposition Principle", "Electric Field E = F/q", "Field due to dipole (axial & equatorial)", "Electric Flux Φ = E A cosθ", "Gauss's Law: ΦE = q_enc/ε₀", "Field due to infinite line charge, sheet, sphere"],
                "formulas": ["F = kq₁q₂/r²", "E = kq/r²", "E_dipole(axial) = 2kp/r³", "E_dipole(equatorial) = kp/r³", "Φ = q/ε₀ (Gauss's Law)"],
                "applications": ["Lightning rods", "Xerography (photocopying)", "Electrostatic precipitators"],
                "common_mistakes": ["Forgetting Gaussian surface must be closed for Gauss's Law", "Direction errors in dipole field"],
                "importance": "15-18% JEE weightage. Gauss's Law problems and dipole questions are JEE staples."
            },
            "Electrostatic Potential & Capacitance": {
                "difficulty": "Hard",
                "description": "Electrostatic potential is the work done per unit charge in bringing a test charge from infinity to a given point. It is scalar, making calculations simpler. Capacitors store energy in electric fields. This chapter covers potential due to charges and dipoles, equipotential surfaces, dielectrics, and combinations of capacitors — all heavily tested in JEE.",
                "roadmap": ["Electric Potential & Relation to E", "Potential due to a Dipole", "Equipotential Surfaces", "Capacitors & Capacitance", "Dielectrics & Energy Stored"],
                "topics": ["V = kq/r (point charge)", "V due to dipole: V = kp cosθ/r²", "Relation: E = −dV/dr", "Capacitance C = Q/V", "Parallel plate: C = ε₀A/d", "Series: 1/C = 1/C₁ + 1/C₂", "Parallel: C = C₁ + C₂", "Energy U = ½CV²", "Dielectric constant K"],
                "formulas": ["V = kq/r", "E = −∇V", "C = ε₀A/d", "U = ½CV²", "U = ε₀E²/2 (energy density)"],
                "applications": ["Camera flash (capacitor discharge)", "Power factor correction", "Touchscreens"],
                "common_mistakes": ["Confusing potential (scalar) with field (vector)", "Mixing up series/parallel capacitor formulas"],
                "importance": "High JEE weightage. Capacitor combinations and energy calculations appear every year."
            },
            "Current Electricity": {
                "difficulty": "Medium",
                "description": "Current Electricity deals with the physics of electric current. Ohm's Law, resistance, resistivity, and temperature dependence form the core. Kirchhoff's Laws (KVL and KCL) are tools for solving complex networks. Wheatstone bridge, Meter Bridge, and Potentiometer appear in every JEE paper.",
                "roadmap": ["Electric Current & Drift Velocity", "Ohm's Law & Resistance", "Temperature Dependence", "Kirchhoff's Laws", "Wheatstone Bridge & Potentiometer"],
                "topics": ["Drift velocity vd = I/nAe", "Resistivity ρ, R = ρl/A", "KVL: ΣV = 0 in loop", "KCL: ΣI = 0 at node", "Wheatstone condition: P/Q = R/S", "EMF & internal resistance: V = E − Ir", "Maximum power transfer: r = R_external"],
                "formulas": ["I = nAevd", "V = IR", "R = ρl/A", "Rₜ = R₀(1 + αT)", "P = VI = I²R = V²/R"],
                "applications": ["Home wiring circuits", "Battery charging systems", "Resistance thermometers"],
                "common_mistakes": ["Wrong sign in KVL loop equations", "Forgetting internal resistance reduces terminal voltage"],
                "importance": "15% JEE weightage. Kirchhoff's Laws and Wheatstone bridge are annual JEE questions."
            },
            "Moving Charges & Magnetism": {
                "difficulty": "Hard",
                "description": "This chapter explores how moving charges create magnetic fields and how magnetic fields exert forces on moving charges. Biot-Savart Law gives field due to a current element; Ampere's Law provides a simpler method for symmetric cases. The Lorentz force causes circular/helical motion — the basis of cyclotrons and mass spectrometers.",
                "roadmap": ["Biot-Savart Law", "Ampere's Circuital Law", "Force on a Charge & Current in B", "Torque on Current Loop", "Moving Coil Galvanometer"],
                "topics": ["Field due to straight wire: B = μ₀I/2πr", "Field at center of circular loop: B = μ₀I/2R", "Solenoid: B = μ₀nI", "Lorentz force: F = q(v×B)", "Cyclotron frequency: f = qB/2πm", "Torque: τ = NIAB sinθ"],
                "formulas": ["B = μ₀I/2πr (wire)", "B = μ₀nI (solenoid)", "F = qvB sinθ", "F = BIl sinθ", "r = mv/qB (circular motion)"],
                "applications": ["Electric motors", "Cyclotron particle accelerator", "MRI machines"],
                "common_mistakes": ["Confusing direction of magnetic force (use Fleming's left-hand rule)", "Forgetting sin θ factor"],
                "importance": "15% JEE weightage. Field calculations and force on conductors are common."
            },
            "Magnetism & Matter": {
                "difficulty": "Medium",
                "description": "Magnetism & Matter studies the magnetic properties of materials — ferromagnetics, paramagnetics, and diamagnetics. Earth's magnetic field elements (declination, dip, horizontal component) and the hysteresis loop of ferromagnetic materials are also covered.",
                "roadmap": ["Bar Magnet as Magnetic Dipole", "Earth's Magnetism", "Magnetic Properties of Materials", "Hysteresis Loop", "Permanent Magnets & Electromagnets"],
                "topics": ["Magnetic moment m = IA", "Field of bar magnet (axial & equatorial)", "Earth's magnetic elements", "Dia, Para, Ferromagnetism", "Curie Temperature", "B-H curve & Hysteresis loss", "Retentivity & Coercivity"],
                "formulas": ["B_axial = μ₀/4π × 2M/r³", "B_equatorial = μ₀/4π × M/r³", "U = −mB cosθ", "χ_m = M/H"],
                "applications": ["Permanent magnets in speakers", "Transformers (low hysteresis loss)", "Magnetic recording media"],
                "common_mistakes": ["Confusing B with H", "Mixing up dia/para/ferromagnetic properties"],
                "importance": "5-8% JEE weightage. Earth's magnetism and classification appear in JEE Main."
            },
            "Electromagnetic Induction": {
                "difficulty": "Hard",
                "description": "Faraday's Law states that changing magnetic flux induces an EMF. Lenz's Law gives its direction. Motional EMF, self-inductance, mutual inductance, and energy in an inductor are key. AC generators, transformers, and eddy currents are tested every year in JEE.",
                "roadmap": ["Magnetic Flux & Faraday's Law", "Lenz's Law", "Motional EMF", "Self & Mutual Inductance", "AC Generator & Transformer"],
                "topics": ["Φ = BA cosθ", "ε = −dΦ/dt (Faraday)", "Lenz's Law — direction of induced current", "Motional EMF ε = Blv", "Self inductance: ε = −L(dI/dt)", "Energy in inductor: U = ½LI²", "Transformer: Vs/Vp = Ns/Np"],
                "formulas": ["ε = −NdΦ/dt", "ε = Blv", "L_solenoid = μ₀n²Al", "U = ½LI²", "Vs/Vp = Ns/Np = Ip/Is"],
                "applications": ["Electric generators & alternators", "Transformers in power grid", "Induction cooktops (eddy currents)"],
                "common_mistakes": ["Sign error in Lenz's law direction", "Confusing L (self) with M (mutual) inductance"],
                "importance": "18-20% JEE weightage. Highest scoring chapter in Class 12 Physics for JEE."
            },
            "Alternating Current": {
                "difficulty": "Hard",
                "description": "AC circuits introduce RMS values, impedance of R, L, C elements, phasors, resonance in LCR circuits, and power factor. The LCR series circuit at resonance is particularly important and tested every year.",
                "roadmap": ["AC Voltage & Current — RMS Values", "Phasors & Phase Relationships", "Reactance of L and C", "LCR Series Circuit & Resonance", "Power in AC Circuits"],
                "topics": ["V_rms = V₀/√2", "Reactance: X_L = ωL, X_C = 1/ωC", "Impedance Z = √(R²+(XL−XC)²)", "Resonance: ω₀ = 1/√LC, Z = R (minimum)", "Quality factor Q = ω₀L/R", "Power factor cos φ = R/Z"],
                "formulas": ["Z = √(R² + (XL−XC)²)", "ω₀ = 1/√(LC)", "Q = ω₀L/R", "P_avg = V_rms I_rms cosφ", "V_rms = V₀/√2"],
                "applications": ["Radio tuning (resonance)", "Power factor correction", "AC power transmission"],
                "common_mistakes": ["Confusing peak and RMS values", "Forgetting resonance condition XL = XC"],
                "importance": "12-15% JEE weightage. Resonance, impedance, and power factor appear every year."
            },
            "Electromagnetic Waves": {
                "difficulty": "Easy",
                "description": "EM Waves are transverse waves with oscillating E and B fields perpendicular to each other and to propagation. Maxwell's equations predicted their existence. The EM spectrum spans from radio waves to gamma rays. Each type's sources and uses are commonly tested in JEE.",
                "roadmap": ["Maxwell's Equations & Displacement Current", "Properties of EM Waves", "Electromagnetic Spectrum", "Speed of EM Waves"],
                "topics": ["Displacement current I_d = ε₀ dΦE/dt", "E & B in phase, perpendicular to propagation", "c = E₀/B₀ = 1/√(μ₀ε₀)", "γ-rays, X-rays, UV, Visible, IR, Microwaves, Radio waves"],
                "formulas": ["c = E/B = 1/√(μ₀ε₀) ≈ 3×10⁸ m/s", "u = ½ε₀E² + B²/2μ₀", "I_d = ε₀ dΦE/dt"],
                "applications": ["Radio/TV broadcasting", "X-ray medical imaging", "Microwave cooking"],
                "common_mistakes": ["Confusing frequency ordering in EM spectrum", "Forgetting E and B are in phase in EM waves"],
                "importance": "5-8% JEE weightage. Know the EM spectrum order."
            },
            "Ray Optics & Optical Instruments": {
                "difficulty": "Medium",
                "description": "Ray Optics treats light as rays and studies reflection and refraction through mirrors and lenses. Total Internal Reflection is the basis of optical fibres. Optical instruments — microscopes and telescopes — are built from lens/mirror combinations. Prisms and dispersion are also important.",
                "roadmap": ["Reflection — Mirrors & Image Formation", "Refraction — Snell's Law & TIR", "Prism & Dispersion", "Lenses", "Optical Instruments"],
                "topics": ["Mirror formula: 1/f = 1/v + 1/u", "Magnification m = −v/u", "Snell's Law", "TIR: sin C = 1/μ", "Lens formula: 1/f = 1/v − 1/u", "Lens maker's equation", "Power P = 1/f", "Microscope & telescope magnification"],
                "formulas": ["1/v − 1/u = 1/f (lens)", "1/v + 1/u = 1/f (mirror)", "μ = c/v = sin i/sin r", "P = P₁ + P₂ (combined lenses)"],
                "applications": ["Optical fibres (TIR)", "Cameras & projectors", "Corrective lenses (spectacles)"],
                "common_mistakes": ["Mixing sign conventions for mirrors vs lenses", "Using mirror formula for lens"],
                "importance": "12-15% JEE weightage. Lens/mirror combinations and TIR are very common."
            },
            "Wave Optics": {
                "difficulty": "Hard",
                "description": "Wave Optics explains interference (YDSE), diffraction, and polarisation. Young's Double Slit Experiment is one of the most important experiments in physics. Brewster's angle for polarisation by reflection and Malus's Law are key concepts.",
                "roadmap": ["Huygens' Principle & Wavefronts", "Young's Double Slit — Interference", "Single Slit Diffraction", "Resolving Power", "Polarisation"],
                "topics": ["Coherent sources for sustained interference", "Path difference Δ = yd/D", "Fringe width β = λD/d", "Constructive: Δ = nλ", "Single slit minima: a sinθ = nλ", "Brewster's angle: tan θ_B = μ", "Malus's Law: I = I₀ cos²θ"],
                "formulas": ["β = λD/d", "y_n = nλD/d", "tan θ_B = μ (Brewster)", "I = I₀ cos²θ (Malus)"],
                "applications": ["Anti-reflection coating", "Holography", "Polaroid sunglasses"],
                "common_mistakes": ["Confusing path difference with phase difference", "Forgetting fringe width is independent of order n"],
                "importance": "12-15% JEE weightage. YDSE is the most tested optics topic."
            },
            "Dual Nature of Radiation & Matter": {
                "difficulty": "Hard",
                "description": "Light behaves as both wave and particle (photon), and matter (electrons) can behave as waves. Einstein's Photoelectric Effect showed light comes in quanta hν. de Broglie hypothesised moving particles have wavelength λ = h/p, confirmed by electron diffraction.",
                "roadmap": ["Photoelectric Effect — Einstein's Theory", "Photon & its Properties", "Threshold Frequency & Work Function", "de Broglie Hypothesis", "Davisson-Germer Experiment"],
                "topics": ["Photon energy E = hν = hc/λ", "Photoelectric equation: KE_max = hν − φ", "Threshold frequency ν₀ = φ/h", "Stopping potential eV₀ = KE_max", "de Broglie wavelength λ = h/p = h/mv"],
                "formulas": ["E = hν", "KE_max = hν − φ", "λ = h/mv (de Broglie)", "λ = h/√(2meV) (electron in field V)"],
                "applications": ["Solar cells", "Photomultiplier tubes", "Electron microscopes"],
                "common_mistakes": ["Confusing frequency with wavelength for threshold", "Forgetting stopping potential ≠ KE directly"],
                "importance": "10-12% JEE weightage. Photoelectric effect and de Broglie wavelength appear every year."
            },
            "Atoms": {
                "difficulty": "Medium",
                "description": "The Atoms chapter traces atomic models: Thomson's plum pudding, Rutherford's nuclear model, and Bohr's quantised model. Bohr's model explains hydrogen line spectra through energy level transitions — Lyman, Balmer, Paschen series.",
                "roadmap": ["Thomson's & Rutherford's Models", "Bohr's Postulates & Energy Levels", "Hydrogen Spectrum — Spectral Series", "Limitations of Bohr's Model"],
                "topics": ["Rutherford's nuclear model", "Bohr's postulates: mvr = nh/2π", "Energy of nth orbit: Eₙ = −13.6/n² eV", "Radius: rₙ = n²a₀", "Spectral series: Lyman(UV), Balmer(visible), Paschen(IR)"],
                "formulas": ["Eₙ = −13.6/n² eV", "rₙ = 0.529 n² Å", "1/λ = R_H(1/n₁² − 1/n₂²)", "v = 2.18×10⁶/n m/s"],
                "applications": ["Atomic clocks", "Laser (stimulated emission)", "Hydrogen fuel cells"],
                "common_mistakes": ["Confusing emission with absorption", "Wrong series — Lyman ends at n=1, Balmer at n=2"],
                "importance": "8-10% JEE weightage. Spectral series and Bohr energy calculations are frequent questions."
            },
            "Nuclei": {
                "difficulty": "Hard",
                "description": "The Nuclei chapter covers nuclear composition, mass defect, binding energy, radioactive decay (α, β, γ), the exponential decay law, half-life, and nuclear fission/fusion — the basis of nuclear reactors and the Sun.",
                "roadmap": ["Nucleus — Size, Mass, Composition", "Mass Defect & Binding Energy", "Radioactivity — α, β, γ Decay", "Radioactive Decay Law & Half-Life", "Nuclear Fission & Fusion"],
                "topics": ["Mass defect Δm, BE = Δmc²", "BE per nucleon curve", "Activity A = λN", "Half life T½ = 0.693/λ", "N = N₀(½)^(t/T½)", "Conservation laws in decay", "Chain reaction & critical mass"],
                "formulas": ["BE = (Zm_p + Nm_n − M)c²", "N = N₀ e^(−λt)", "T½ = 0.693/λ", "A = λN", "Q = (m_reactants − m_products)c²"],
                "applications": ["Nuclear power plants (fission)", "Carbon dating (radioactive decay)", "Hydrogen bomb (fusion)"],
                "common_mistakes": ["Confusing half-life with mean life (τ = T½/0.693)", "Forgetting to use 931.5 MeV/amu"],
                "importance": "10-12% JEE weightage. Binding energy curve and decay calculations appear every year."
            },
            "Semiconductor Electronics": {
                "difficulty": "Medium",
                "description": "Semiconductor Electronics covers energy bands, doping (p-type/n-type), p-n junction diodes, rectifiers, Zener voltage regulator, transistors (BJT in CE), and logic gates — the foundation of all digital electronics.",
                "roadmap": ["Energy Bands — Conductor, Insulator, Semiconductor", "p-type & n-type Semiconductors", "p-n Junction Diode & Rectification", "Transistor — CE Configuration", "Logic Gates & Boolean Algebra"],
                "topics": ["Intrinsic & extrinsic semiconductors", "Depletion layer & barrier potential", "Half-wave and full-wave rectifier", "Zener diode as voltage regulator", "BJT: β = Ic/Ib", "NAND & NOR as universal gates", "Truth tables"],
                "formulas": ["I = I₀(e^(eV/kT) − 1) (diode)", "β = Ic/Ib", "α = Ic/Ie", "β = α/(1−α)"],
                "applications": ["Mobile phones & computers", "Solar cells (p-n junction)", "LED lighting"],
                "common_mistakes": ["Confusing α and β for transistor", "Forgetting NAND = NOT AND (outputs inverted)"],
                "importance": "10% JEE weightage. Logic gates and diode rectification are easy marks in JEE Main."
            },
            "Communication Systems": {
                "difficulty": "Easy",
                "description": "Communication Systems covers how information is transmitted using EM signals. Modulation (AM/FM), demodulation, bandwidth, signal propagation modes, and internet/mobile/satellite communication are key topics — mostly conceptual, easy scoring.",
                "roadmap": ["Elements of Communication System", "Signal Bandwidth & Channel", "Modulation — AM & FM", "Signal Propagation Modes", "Internet, Mobile & Satellite Communication"],
                "topics": ["Source, Transmitter, Channel, Receiver, Sink", "Bandwidth of AM signal = 2f_m", "Modulation index m = Am/Ac", "Ground wave, Sky wave, Space wave propagation", "Maximum usable frequency (MUF)"],
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
                "topics": ["Galvanic cell: anode (oxidation, -) and cathode (reduction, +)", "Salt bridge: maintains electrical neutrality", "Standard hydrogen electrode (SHE): E° = 0 V", "E°cell = E°cathode - E°anode", "Electrochemical series — relative reactivity", "Nernst equation: E = E° - (0.0591/n)log Q at 25°C", "ΔG = -nFE; ΔG° = -nFE° = -RT ln K", "Faraday's 1st law: m = ZIt = (M/nF)It", "Faraday's 2nd law: m₁/m₂ = E₁/E₂", "Kohlrausch's Law: Λ°m = Σλ°(cations) + Σλ°(anions)", "Dry cell: Zn-MnO₂", "Lead storage battery: Pb-PbO₂-H₂SO₄", "Fuel cells: H₂-O₂ cell"],
                "formulas": ["E°cell = E°cathode - E°anode", "E = E° - (0.0591/n)log Q (Nernst)", "ΔG = -nFE", "m = (M × I × t)/(n × F) (Faraday)", "Λm = κ × 1000/c"],
                "applications": ["Lead-acid batteries in vehicles", "Electroplating (silver, gold)", "Fuel cells in spacecraft"],
                "common_mistakes": ["Confusing anode/cathode sign in different cell types", "Forgetting temperature in Nernst equation", "Wrong formula for Faraday's Law"],
                "importance": "15-18% JEE weightage. Nernst equation and electrolysis problems appear every year."
            },
            "Chemical Kinetics": {
                "difficulty": "Hard",
                "description": "Chemical Kinetics studies the rate (speed) of chemical reactions and the factors that affect it. The rate law relates reaction rate to reactant concentrations. Integrated rate laws allow calculation of concentration at any time. Arrhenius Equation connects rate constant to temperature and activation energy. Catalysts speed up reactions by providing an alternative pathway with lower activation energy.",
                "roadmap": ["Rate of Reaction & Factors Affecting It", "Rate Law & Order of Reaction", "Integrated Rate Equations", "Temperature Dependence — Arrhenius", "Collision Theory & Catalysis"],
                "topics": ["Rate law: Rate = k[A]^m[B]^n", "Order = m + n (from experiment, not stoichiometry)", "Zero order: t½ = [A]₀/2k", "First order: ln[A] = ln[A]₀ - kt; t½ = 0.693/k", "Units of k: zero order (mol L⁻¹ s⁻¹), first order (s⁻¹)", "Arrhenius: k = Ae^(-Ea/RT)", "log(k₂/k₁) = Ea/2.303R × (T₂-T₁)/(T₁T₂)", "Pseudo-first order reactions"],
                "formulas": ["Rate = k[A]^m[B]^n", "t½ = 0.693/k (first order)", "k = Ae^(-Ea/RT) (Arrhenius)", "log(k₂/k₁) = Ea/2.303R × (T₂-T₁)/(T₁T₂)"],
                "applications": ["Refrigeration (slowing food spoilage)", "Industrial catalyst design", "Drug metabolism rates in pharmacology"],
                "common_mistakes": ["Determining order from balanced equation (order is from experiment!)", "Wrong units of k for different orders", "Confusing t½ formulas for zero, first, second order"],
                "importance": "12-15% JEE weightage. Graph-based questions (ln[A] vs t, [A] vs t) are very common."
            },
            "Surface Chemistry": {
                "difficulty": "Medium",
                "description": "Surface Chemistry deals with phenomena occurring at surfaces and interfaces. Adsorption explains how substances accumulate at surfaces (physosorption vs chemisorption). Colloids are heterogeneous mixtures with particle sizes 1-1000 nm that show unique properties like the Tyndall effect and Brownian motion.",
                "roadmap": ["Adsorption — Physosorption vs Chemisorption", "Freundlich & Langmuir Adsorption Isotherms", "Catalysis", "Colloidal State — Types & Properties", "Emulsions & Applications"],
                "topics": ["Physosorption: weak van der Waals, reversible, multilayer", "Chemisorption: strong chemical bond, irreversible, monolayer", "Freundlich adsorption isotherm: x/m = kp^(1/n)", "Tyndall effect: scattering of light by colloidal particles", "Brownian motion: random movement due to molecular bombardment", "Electrophoresis: movement of colloidal particles in electric field", "Hardy-Schulze rule: higher valency = faster coagulation", "Types of colloids: lyophilic (stable) and lyophobic (unstable)", "Micelles: soap/detergent in water above CMC"],
                "formulas": ["x/m = kp^(1/n) (Freundlich)", "x/m = ap/(1+bp) (Langmuir)"],
                "applications": ["Catalytic converters (heterogeneous catalysis)", "Dialysis for kidney patients", "Food emulsions (mayonnaise, milk)"],
                "common_mistakes": ["Confusing lyophilic and lyophobic colloids", "Forgetting Tyndall effect is shown by colloids NOT true solutions"],
                "importance": "8-10% JEE weightage. Mostly conceptual — 2-3 one-liner questions in JEE Main."
            },
            "d & f Block Elements": {
                "difficulty": "Hard",
                "description": "d-Block elements (transition metals, Groups 3-12) have partially filled d-orbitals giving them unique properties: variable oxidation states, colored compounds, magnetic properties, and catalytic activity. f-Block elements (lanthanoids and actinoids) have 4f and 5f orbitals being filled. K₂Cr₂O₇ and KMnO₄ are the most important oxidizing agents tested in JEE.",
                "roadmap": ["Electronic Configuration — (n-1)d¹⁻¹⁰ns¹⁻²", "Properties of Transition Metals", "Important Compounds: K₂Cr₂O₇ & KMnO₄", "Lanthanoids & Lanthanoid Contraction", "Actinoids & Comparison"],
                "topics": ["Variable valency: multiple d-electrons can be involved in bonding", "Color due to d-d transitions in unpaired electrons", "Magnetic moment: μ = √n(n+2) BM (n = unpaired e⁻)", "KMnO₄: Mn goes from +7 to +2 (acid), +4 (neutral), +6 (base)", "K₂Cr₂O₇: Cr goes from +6 to +3 in acidic reactions", "Lanthanoid contraction: steady decrease in atomic radius due to poor 4f shielding", "Consequence: Zr ≈ Hf (same radius), similar chemistry"],
                "formulas": ["μ = √n(n+2) BM (magnetic moment)"],
                "applications": ["Fe catalyst in Haber process", "Pt/Pd in catalytic converters", "Ti alloys in aerospace"],
                "common_mistakes": ["Wrong electronic config of Cr (3d⁵4s¹) and Cu (3d¹⁰4s¹)", "Forgetting KMnO₄ product changes with medium (acidic/neutral/basic)"],
                "importance": "12-15% JEE weightage. KMnO₄/K₂Cr₂O₇ reactions and magnetic moment appear every year."
            },
            "Coordination Compounds": {
                "difficulty": "Hard",
                "description": "Coordination Compounds are compounds where a central metal atom/ion is bonded to surrounding ligands via coordinate (dative) bonds. Werner's Theory first explained their structure. IUPAC nomenclature, types of isomerism (geometrical, optical, ionization, linkage), Valence Bond Theory (VBT) for hybridization, and Crystal Field Theory (CFT) for color and magnetism are all tested extensively.",
                "roadmap": ["Werner's Theory & Terminology", "IUPAC Nomenclature", "Isomerism in Coordination Compounds", "Valence Bond Theory (VBT)", "Crystal Field Theory (CFT)"],
                "topics": ["Denticity: monodentate, bidentate, polydentate, chelate", "IUPAC: ligands (alphabetical), then metal, then ox. state", "Geometrical isomerism: cis-trans in square planar and octahedral", "Optical isomerism: non-superimposable mirror images (chirality)", "Ionization isomerism: [Co(NH₃)₅Br]SO₄ vs [Co(NH₃)₅SO₄]Br", "VBT: hybridization determines geometry", "Strong field ligands (CO, CN⁻): large Δ, low spin", "Weak field ligands (F⁻, Cl⁻): small Δ, high spin"],
                "formulas": ["EAN = Z - oxidation state + 2 × CN (18-electron rule)"],
                "applications": ["Haemoglobin: Fe²⁺ complex with O₂", "Cisplatin [Pt(NH₃)₂Cl₂]: cancer drug", "EDTA: chelation therapy"],
                "common_mistakes": ["Wrong IUPAC naming order (ligands first alphabetically, then metal)", "Confusing geometrical and optical isomerism conditions"],
                "importance": "12-15% JEE weightage. IUPAC naming, isomerism, and CFT color are annual JEE topics."
            },
            "Haloalkanes & Haloarenes": {
                "difficulty": "Hard",
                "description": "Haloalkanes are alkyl halides (RX) formed by replacing H with a halogen. Haloarenes are aryl halides (ArX). SN1 and SN2 are competing nucleophilic substitution mechanisms. Key named reactions (Finkelstein, Swarts, Wurtz-Fittig, Grignard reagent formation) are frequently tested in JEE.",
                "roadmap": ["Classification & Nomenclature", "Preparation Methods", "SN1 vs SN2 Mechanisms", "Elimination Reactions (E1/E2)", "Reactions & Named Reactions"],
                "topics": ["SN2: bimolecular, backside attack, inversion (Walden inversion), 1° halides prefer", "SN1: unimolecular, carbocation intermediate, racemization, 3° halides prefer", "Reactivity: RI > RBr > RCl > RF", "Grignard reagent: RX + Mg → RMgX (dry ether)", "Finkelstein reaction: RI exchange with NaI in acetone", "Swarts reaction: RF from RCl + AgF or SbF₃", "Wurtz reaction: RX + 2Na → R-R", "Haloarenes: less reactive toward nucleophilic substitution"],
                "formulas": ["SN2: rate = k[RX][Nu⁻] (bimolecular)", "SN1: rate = k[RX] (unimolecular)"],
                "applications": ["Chloroform (CHCl₃) — anesthetic", "DDT — pesticide (banned)", "Freons (CFCs) — refrigerants (banned — ozone depletion)"],
                "common_mistakes": ["Wrong mechanism for 3° vs 1° substrates", "Forgetting Grignard reagent requires anhydrous (no water)"],
                "importance": "12-15% JEE weightage. SN1/SN2 mechanism is JEE staple — appears every year."
            },
            "Alcohols, Phenols & Ethers": {
                "difficulty": "Hard",
                "description": "Alcohols (ROH) and phenols (ArOH) are hydroxyl compounds with contrasting properties — phenols are much more acidic than alcohols due to resonance stabilization of phenoxide ion. Lucas test distinguishes 1°, 2°, 3° alcohols; Reimer-Tiemann gives o-hydroxybenzaldehyde from phenol.",
                "roadmap": ["Classification & Nomenclature", "Preparation of Alcohols & Phenols", "Physical Properties — H-bonding", "Chemical Reactions", "Phenol Reactions & Acidity"],
                "topics": ["Lucas test: ZnCl₂/HCl — 3° (instant), 2° (5 min), 1° (no reaction)", "Dehydration: 1° → E2 (180°C), 2°/3° → E1 (low temp)", "Victor Meyer test: 1° → red, 2° → blue, 3° → colorless", "Phenol acidity: pKa ~10 — more acidic than alcohol", "EAS at ortho and para positions (+M effect of OH)", "Reimer-Tiemann: CHCl₃ + NaOH + Phenol → o-hydroxybenzaldehyde", "Kolbe's reaction: Na-phenoxide + CO₂ → salicylic acid", "Williamson synthesis: R-O⁻ + R'X → R-O-R' (SN2)"],
                "formulas": ["Williamson synthesis: R-O⁻ + R'X → R-O-R' (SN2)"],
                "applications": ["Ethanol as biofuel", "Phenol in disinfectants (Dettol)", "Aspirin from salicylic acid"],
                "common_mistakes": ["Forgetting Lucas test only works for 3° immediately", "Wrong product for Victor Meyer test", "Confusing Reimer-Tiemann with Kolbe reaction"],
                "importance": "15-18% JEE weightage. One of the highest scoring organic chapters."
            },
            "Aldehydes, Ketones & Carboxylic Acids": {
                "difficulty": "Hard",
                "description": "Carbonyl compounds (C=O) are among the most reactive functional groups. Aldehydes (RCHO) are more reactive than ketones (RCOR') toward nucleophilic addition. Named reactions — Aldol, Cannizzaro, Clemmensen, Wolff-Kishner, Tollens', Fehling's, Haloform, HVZ — are extremely important for JEE.",
                "roadmap": ["Carbonyl Group — Structure & Reactivity", "Nucleophilic Addition Reactions", "Oxidation & Reduction", "Named Reactions (Aldol, Cannizzaro, Clemmensen)", "Carboxylic Acids & Derivatives"],
                "topics": ["Nucleophilic addition: HCN, RMgX (Grignard), NaHSO₃, alcohols → hemiacetals/acetals", "Reactivity: HCHO > RCHO > RCOR'", "Aldol condensation: α-H + NaOH → β-hydroxy carbonyl compound", "Cannizzaro reaction: HCHO + NaOH (no α-H) → alcohol + acid", "Clemmensen reduction: Zn/Hg + HCl → alkane (acidic)", "Wolff-Kishner reduction: H₂NNH₂ + KOH → alkane (basic)", "Haloform reaction: RCOCH₃ + X₂ + NaOH → CHX₃", "Tollens' test: all aldehydes; Fehling's: aliphatic aldehydes only"],
                "formulas": ["Aldol: 2CH₃CHO → CH₃CH(OH)CH₂CHO"],
                "applications": ["Formaldehyde (HCHO) in formalin (preservative)", "Acetic acid in vinegar", "Aspirin from salicylic acid esterification"],
                "common_mistakes": ["Cannizzaro requires no α-H — else Aldol occurs", "Tollens' test works for all aldehydes; Fehling's only for aliphatic", "Confusing Clemmensen (acid) with Wolff-Kishner (base)"],
                "importance": "18-20% JEE weightage. Highest scoring organic chapter — named reactions are mandatory."
            },
            "Amines": {
                "difficulty": "Hard",
                "description": "Amines are nitrogen-containing organic bases. Classified as primary (1°), secondary (2°), tertiary (3°). Basicity: aliphatic > ammonia > aniline. Diazonium salts formed from primary aryl amines are key synthetic intermediates for azo dyes, Sandmeyer reactions, and coupling reactions.",
                "roadmap": ["Classification & Nomenclature", "Preparation Methods", "Basicity of Amines", "Reactions", "Diazonium Salts & Coupling"],
                "topics": ["Preparation: reduction of nitro compounds, Gabriel phthalimide (1° amines only), Hoffmann bromamide (reduces C by 1)", "Basicity: 2° aliphatic > 1° aliphatic > 3° aliphatic > NH₃ > aniline", "Hinsberg's test: distinguish 1°, 2°, 3° amines", "Carbylamine test: 1° amines + CHCl₃ + NaOH → isocyanide", "Diazonium salts: ArNH₂ + NaNO₂ + HCl (0-5°C) → ArN₂⁺Cl⁻", "Sandmeyer: ArN₂⁺ + CuCN → ArCN, + CuCl → ArCl, + CuBr → ArBr", "Coupling reaction: ArN₂⁺ + ArNR₂ → azo dye (orange-red)"],
                "formulas": ["Basicity: 2°(aliphatic) > 1° > 3° > NH₃ > aniline"],
                "applications": ["Aniline — precursor for dyes, drugs, explosives", "Diazonium coupling — azo dyes in textiles"],
                "common_mistakes": ["Forgetting carbylamine test is only for primary amines", "Diazonium salt preparation requires 0-5°C — higher temp decomposes it"],
                "importance": "12-15% JEE weightage. Named reactions (Hoffmann, Sandmeyer, carbylamine) appear every year."
            },
            "Biomolecules": {
                "difficulty": "Medium",
                "description": "Biomolecules are the large organic molecules that form the basis of all living systems. Carbohydrates provide energy and structural support, proteins are functional molecules (enzymes, antibodies), nucleic acids (DNA, RNA) store and transmit genetic information.",
                "roadmap": ["Carbohydrates — Classification & Glucose", "Proteins — Amino Acids & Peptides", "Nucleic Acids — DNA & RNA", "Vitamins & Enzymes", "Lipids"],
                "topics": ["Monosaccharides: glucose (C₆H₁₂O₆), fructose (C₆H₁₂O₆)", "Disaccharides: sucrose (non-reducing), maltose/lactose (reducing)", "Polysaccharides: starch (food), cellulose (structural), glycogen (animal)", "Reducing sugars: give Fehling's/Tollens' test", "Peptide bond: -CO-NH- linkage", "Primary → Secondary (α-helix, β-sheet) → Tertiary → Quaternary structure", "Denaturation: loss of structure without breaking peptide bonds", "DNA: double helix, A=T (2 H-bonds), G≡C (3 H-bonds)", "RNA: single strand, uracil replaces thymine", "Vitamins: fat-soluble (A, D, E, K) and water-soluble (B, C)"],
                "formulas": ["Peptide bond: -CO-NH-", "n amino acids → n-1 peptide bonds + n-1 H₂O"],
                "applications": ["Insulin (protein hormone) for diabetes", "DNA fingerprinting in forensics"],
                "common_mistakes": ["Sucrose is non-reducing (no free OH on anomeric C)", "Denaturation ≠ hydrolysis of peptide bonds", "RNA has uracil; DNA has thymine"],
                "importance": "8-10% JEE weightage. Glucose structure, DNA base pairing, reducing sugar tests are frequent."
            },
            "Polymers": {
                "difficulty": "Easy",
                "description": "Polymers are giant molecules formed by linking many small repeating units called monomers. Addition polymerization involves alkenes (no byproduct), while condensation polymerization forms a small molecule as byproduct. Natural rubber is a polymer of isoprene; vulcanization adds sulfur crosslinks.",
                "roadmap": ["Classification of Polymers", "Addition Polymerization (Chain Growth)", "Condensation Polymerization (Step Growth)", "Natural & Synthetic Rubber", "Biodegradable Polymers"],
                "topics": ["Addition polymers: polythene, PVC, teflon, polystyrene", "HDPE (linear) vs LDPE (branched)", "Condensation polymers: nylon-6,6, nylon-6, terylene/dacron", "Bakelite: phenol + formaldehyde, thermosetting, cross-linked", "Natural rubber: cis-1,4-polyisoprene", "Vulcanization: S₈ creates crosslinks", "Buna-S: butadiene + styrene (SBR)"],
                "formulas": ["Degree of polymerization DP = MW(polymer)/MW(monomer)"],
                "applications": ["LDPE in plastic bags", "Nylon-6,6 in parachutes and stockings", "Bakelite in electrical insulators"],
                "common_mistakes": ["Confusing nylon-6 (one monomer) with nylon-6,6 (two monomers)", "Terylene is a polyester, not polyamide"],
                "importance": "6-8% JEE weightage. Monomers of common polymers and Bakelite structure are tested."
            },
            "Chemistry in Everyday Life": {
                "difficulty": "Easy",
                "description": "Chemistry in Everyday Life explores practical applications in medicine, food, and hygiene. Covers analgesics, antipyretics, antibiotics, antiseptics, disinfectants, antacids, tranquilizers, food preservatives, and cleansing agents. Mostly factual — easy scoring for JEE Main.",
                "roadmap": ["Drugs & Medicines Classification", "Analgesics, Antipyretics & Antibiotics", "Antiseptics & Disinfectants", "Food Chemicals", "Cleansing Agents — Soaps & Detergents"],
                "topics": ["Analgesics: non-narcotic (aspirin, paracetamol) and narcotic (morphine)", "Antibiotics: narrow spectrum (penicillin) vs broad spectrum (chloramphenicol)", "Antiseptics: safe for skin (0.2% phenol, iodoform)", "Disinfectants: for inanimate objects (1% phenol, bleach)", "Antacids: NaHCO₃, Mg(OH)₂ — reduce stomach acidity", "Tranquilizers: diazepam/Valium, meprobamate", "Saccharin: 550× sweeter than sugar, no calories", "Soap: RCOONa, saponification", "Detergent: works in hard water (no scum)"],
                "formulas": ["Saponification: fat + NaOH → soap + glycerol"],
                "applications": ["Aspirin: pain + fever + anti-clotting", "Penicillin: first antibiotic", "Detergents in hard water"],
                "common_mistakes": ["Antiseptic ≠ disinfectant (antiseptics are safer for living tissue)", "Soap doesn't work in hard water — detergent does"],
                "importance": "5-7% JEE weightage. One-liner theory questions — easy scoring in JEE Main."
            },
            "Solutions": {
                "difficulty": "Hard",
                "description": "Solutions explores the properties of liquid mixtures. Concentration expressions, Henry's Law, Raoult's Law, ideal/non-ideal solutions, and colligative properties (ΔTf, ΔTb, π) form the core. Van't Hoff factor (i) accounts for association and dissociation of solutes.",
                "roadmap": ["Concentration Expressions", "Henry's Law & Raoult's Law", "Ideal & Non-Ideal Solutions", "Colligative Properties (ΔTf, ΔTb, π)", "Van't Hoff Factor"],
                "topics": ["Molarity M = moles of solute/L", "Molality m = moles of solute/kg solvent", "Henry's Law: p = KH × x (gas dissolved in liquid)", "Raoult's Law: p = p° × mole fraction", "Positive deviation: e.g., acetone-water; Negative: e.g., acetone-CHCl₃", "Azeotropes: constant boiling mixtures", "ΔTf = Kf × m × i", "ΔTb = Kb × m × i", "π = iMRT (osmotic pressure)", "Reverse osmosis: apply pressure > osmotic pressure"],
                "formulas": ["ΔTf = Kf × m", "ΔTb = Kb × m", "π = CRT = MRT", "Van't Hoff: i = observed/theoretical"],
                "applications": ["Antifreeze (ethylene glycol): lowers f.p.", "RO water purification plants", "IV saline: isotonic with blood (0.9% NaCl)"],
                "common_mistakes": ["Molarity changes with temperature, molality doesn't", "Forgetting Van't Hoff factor i for electrolytes"],
                "importance": "12-15% JEE weightage. Colligative property numericals (ΔTf, ΔTb, π) appear every year."
            },
            "Solid State": {
                "difficulty": "Hard",
                "description": "Solid State studies the three-dimensional arrangement of atoms in crystalline solids. Unit cells (SCC, BCC, FCC), packing efficiency, void calculations, and crystal defects (Schottky and Frenkel) are heavily tested numerically in JEE.",
                "roadmap": ["Types of Solids — Crystalline vs Amorphous", "Classification of Crystalline Solids", "Unit Cells & Crystal Lattices", "Packing Efficiency & Void Calculation", "Crystal Defects & Properties"],
                "topics": ["SCC: CN = 6, PE = 52.4%, 1 atom/unit cell", "BCC: CN = 8, PE = 68%, 2 atoms/unit cell (r = √3a/4)", "FCC/CCP: CN = 12, PE = 74%, 4 atoms/unit cell (r = a/2√2)", "In FCC: 8 tetrahedral voids, 4 octahedral voids per unit cell", "NaCl structure: FCC with Na in octahedral voids (CN = 6)", "ZnS (zinc blende): FCC with Zn in half tetrahedral voids (CN = 4)", "Schottky defect: cation + anion vacancies → decreases density", "Frenkel defect: smaller ion moves to interstitial site → no density change"],
                "formulas": ["Density = Z × M / (a³ × NA)", "Packing efficiency (FCC) = 74%", "Packing efficiency (BCC) = 68%"],
                "applications": ["Diamond (covalent crystal) hardest solid", "Semiconductor doping (crystal defects)"],
                "common_mistakes": ["Counting atoms in unit cell (corner = 1/8, face = 1/2, edge = 1/4, body = 1)", "Schottky = vacancy defect; Frenkel = interstitial defect"],
                "importance": "10-12% JEE weightage. Unit cell density formula and packing efficiency are annual JEE numericals."
            }
        },

        "Maths": {
            "Calculus": {
                "difficulty": "Hard",
                "description": "Calculus covers continuity & differentiability, differentiation applications, integration techniques, definite integrals, and area under curves.",
                "roadmap": ["Continuity & Differentiability", "Differentiation Applications", "Integration Techniques", "Definite Integrals", "Area Under Curves"],
                "topics": ["Derivative Rules", "Rolle's & Lagrange's Theorem", "Maxima & Minima", "Integration by Parts", "Definite Integral Properties"],
                "formulas": ["d/dx(xⁿ) = nxⁿ⁻¹", "∫x^n dx = x^(n+1)/(n+1)", "∫udv = uv - ∫vdu", "dy/dx = (dy/dt)/(dx/dt)"],
                "importance": "Highest 25-30% board weightage. Core JEE Advanced topic."
            },
            "Vectors": {
                "difficulty": "Medium",
                "description": "Vectors covers vector algebra basics, dot product, cross product, scalar triple product, and applications in geometry.",
                "roadmap": ["Vector Algebra Basics", "Dot Product", "Cross Product", "Scalar Triple Product", "Applications"],
                "topics": ["Direction Cosines & Ratios", "Position Vector", "Projections", "Area of Parallelogram", "Scalar & Vector Triple Products"],
                "formulas": ["a·b = |a||b|cosθ", "a×b = |a||b|sinθ n̂", "[a b c] = a·(b×c)"],
                "importance": "10-12% board weightage. 3D geometry foundation for JEE."
            },
            "3D Geometry": {
                "difficulty": "Hard",
                "description": "3D Geometry covers direction cosines, equation of line and plane, angle between lines/planes, and distance formulas in 3D space.",
                "roadmap": ["Direction Cosines", "Equation of Line", "Equation of Plane", "Angle Between Lines/Planes", "Distance Formulas"],
                "topics": ["Cartesian & Vector Form", "Skew Lines", "Coplanar Lines", "Perpendicular Distance", "Shortest Distance"],
                "formulas": ["(x-x₁)/a = (y-y₁)/b = (z-z₁)/c", "ax + by + cz + d = 0", "d = |ax₁+by₁+cz₁+d|/√(a²+b²+c²)"],
                "importance": "15% board weightage. 3D visualization needed for JEE."
            },
            "Matrices & Determinants": {
                "difficulty": "Hard",
                "description": "Matrices & Determinants covers matrix operations, types of matrices, determinant properties, adjoint & inverse, and system of equations.",
                "roadmap": ["Matrix Operations", "Types of Matrices", "Determinant Properties", "Adjoint & Inverse", "System of Equations"],
                "topics": ["Addition, Multiplication", "Transpose, Symmetric", "Cofactor Expansion", "Cramer's Rule", "Rank of Matrix"],
                "formulas": ["|AB| = |A||B|", "A⁻¹ = adj(A)/|A|", "A(adj A) = |A|I"],
                "importance": "12-15% board weightage. Matrix problems frequent in JEE."
            },
            "Probability": {
                "difficulty": "Hard",
                "description": "Probability covers random experiments, conditional probability, Bayes' Theorem, random variables, and probability distributions.",
                "roadmap": ["Random Experiments", "Conditional Probability", "Bayes' Theorem", "Random Variables", "Probability Distributions"],
                "topics": ["Addition & Multiplication Theorems", "Independent Events", "Binomial Distribution", "Mean & Variance"],
                "formulas": ["P(A∪B) = P(A) + P(B) - P(A∩B)", "P(A|B) = P(A∩B)/P(B)", "P(E) = ΣP(Ei)P(E|Ei)"],
                "importance": "10-12% board weightage. Combinatorics integration in JEE."
            },
            "Differential Equations": {
                "difficulty": "Hard",
                "description": "Differential Equations covers order & degree, variables separable method, homogeneous equations, linear equations, and growth & decay applications.",
                "roadmap": ["Order & Degree", "Variables Separable", "Homogeneous Equations", "Linear Equations", "Applications"],
                "topics": ["First Order DE", "dy/dx + Py = Q", "Integrating Factor", "Formation of DE", "Growth & Decay"],
                "formulas": ["IF = e^∫Pdx", "y·IF = ∫Q·IF dx"],
                "importance": "8-10% board weightage. Application-based JEE problems."
            },
            "Relations & Functions": {
                "difficulty": "Medium",
                "description": "Relations & Functions covers types of relations and functions, composition, inverse functions, and binary operations.",
                "roadmap": ["Types of Relations", "Types of Functions", "Composition", "Inverse Functions", "Binary Operations"],
                "topics": ["Reflexive, Symmetric, Transitive", "One-One, Onto, Bijective", "fog & gof", "Invertible Functions"],
                "formulas": ["(fog)(x) = f(g(x))", "f⁻¹(f(x)) = x"],
                "importance": "8% board weightage. Conceptual clarity for JEE."
            },
            "Linear Programming": {
                "difficulty": "Easy",
                "description": "Linear Programming covers problem formulation, graphical solution, feasible region, and optimization using the corner point method.",
                "roadmap": ["Problem Formulation", "Graphical Solution", "Feasible Region", "Optimization", "Corner Point Method"],
                "topics": ["Objective Function", "Constraints", "Bounded & Unbounded Regions", "Maximum & Minimum Values"],
                "formulas": ["Optimize Z = ax + by subject to constraints"],
                "importance": "6-8% board weightage. Easy graphical marks."
            },
            "Applications of Integrals": {
                "difficulty": "Medium",
                "description": "Applications of Integrals covers area between curves, area of circle/ellipse, and volume of solids of revolution.",
                "roadmap": ["Area Between Curves", "Area of Circle/Ellipse", "Volume of Solids"],
                "topics": ["Area under y = f(x)", "Area between two curves", "Volume by rotation"],
                "formulas": ["Area = ∫[a to b] f(x)dx", "Area = ∫[a to b] [f(x) - g(x)]dx"],
                "importance": "8-10% board weightage. Visualization needed for JEE."
            }
        }
    }
};

const app = {
    currentSubject: 'Physics',
    currentClass: '11',
    activeChapter: null,
    chapters: [],
    isSidebarActive: false,
    cache: {},
    bookmarks: JSON.parse(localStorage.getItem('bookmarks') || '[]'),

    init() {
        this.setupEventListeners();
        this.loadState();
        this.syncSidebarUI();
        this.showLanding();
        lucide.createIcons();
    },

    refreshContent(classNum) {
        if (!classNum) classNum = this.currentClass;
        this.currentClass = classNum;
        localStorage.setItem('userClass', classNum);
        this.syncSidebarUI();
        this.loadSubject(this.currentSubject);
    },

    loadState() {
        const savedClass = localStorage.getItem('userClass');
        const savedSubject = localStorage.getItem('userSubject');
        const savedChapter = localStorage.getItem('lastChapter');
        if (savedClass) this.currentClass = savedClass;
        if (savedSubject) this.currentSubject = savedSubject;
        if (savedChapter) this.activeChapter = savedChapter;
    },

    syncSidebarUI() {
        document.querySelectorAll('[data-class]').forEach(btn => {
            btn.classList.toggle('active', btn.dataset.class === this.currentClass);
        });
        document.querySelectorAll('[data-subj-btn]').forEach(btn => {
            btn.classList.toggle('active', btn.dataset.subjBtn === this.currentSubject);
        });
    },

    selectSidebarClass(event, classNum) {
        this.currentClass = classNum;
        localStorage.setItem('userClass', classNum);
        this.syncSidebarUI();
        this.loadSubject(this.currentSubject);
    },

    setupEventListeners() {
        window.addEventListener('scroll', () => {
            this.updateProgressBar();
            this.toggleBackToTop();
        });
        window.addEventListener('load', () => {
            const overlay = document.getElementById('entranceOverlay');
            if (overlay) {
                setTimeout(() => overlay.style.opacity = '0', 500);
                setTimeout(() => overlay.style.display = 'none', 1000);
            }
        });
    },

    showLanding() {
        document.getElementById('landingPage').style.display = 'flex';
        document.getElementById('appContainer').style.display = 'none';
        const aiBtn = document.getElementById('shishimanuBtn');
        if (aiBtn) aiBtn.style.display = 'none';
        const aiPanel = document.getElementById('shishimanuPanel');
        if (aiPanel) aiPanel.classList.add('hidden');
    },

    showApp() {
        document.getElementById('landingPage').style.display = 'none';
        document.getElementById('appContainer').style.display = 'flex';
        document.body.style.overflow = 'auto';
        const aiBtn = document.getElementById('shishimanuBtn');
        if (aiBtn) aiBtn.style.display = '';
    },

    showOnboarding() {
        document.getElementById('onboardingModal').classList.remove('hidden');
    },

    selectClass(event, classNum) {
        this.currentClass = classNum;
        document.querySelectorAll('[data-class]').forEach(btn => btn.classList.remove('selected'));
        event.target.closest('button').classList.add('selected');
        this.checkOnboardingComplete();
    },

    selectOnboardingSubject(event, subject) {
        this.currentSubject = subject;
        document.querySelectorAll('[data-subject]').forEach(btn => btn.classList.remove('selected'));
        event.target.closest('button').classList.add('selected');
        this.checkOnboardingComplete();
    },

    checkOnboardingComplete() {
        const btn = document.getElementById('continueBtn');
        btn.disabled = !(this.currentClass && this.currentSubject);
    },

    completeOnboarding() {
        localStorage.setItem('hasOnboarded', 'true');
        localStorage.setItem('userClass', this.currentClass);
        localStorage.setItem('userSubject', this.currentSubject);
        document.getElementById('onboardingModal').classList.add('hidden');
        this.showApp();
        this.loadSubject(this.currentSubject);
    },

    async loadSubject(subject) {
        this.currentSubject = subject;
        localStorage.setItem('userSubject', subject);
        const headerSubjectEl = document.getElementById('headerSubject');
        if (headerSubjectEl) headerSubjectEl.innerText = subject;
        document.querySelectorAll('[data-subj-btn]').forEach(btn => {
            btn.classList.toggle('active', btn.dataset.subjBtn === subject);
        });

        const data = await this.fetchData(this.currentClass, subject);
        if (data && Object.keys(data).length) {
            this.chapters = Object.keys(data);
            this.renderSidebar(data);
            this.renderDashboard(data);
        } else {
            this.chapters = [];
            this.renderSidebar({});
            document.getElementById('chapterContent').innerHTML =
                '<div class="error-msg">No chapters found for this subject/class.</div>';
        }
    },

    async fetchData(classId, subject) {
        const key = `${classId}-${subject}`;
        if (this.cache[key]) return this.cache[key];
        const data = (CHAPTER_DATA[classId] || {})[subject] || {};
        this.cache[key] = data;
        return data;
    },

    renderSidebar(data) {
        const list = document.getElementById('chapterList');
        list.innerHTML = '';
        this.chapters.forEach((name, idx) => {
            const isBookmarked = this.bookmarks.includes(`${this.currentClass}-${this.currentSubject}-${name}`);
            const item = document.createElement('div');
            item.className = `chapter-item ${this.activeChapter === name ? 'active' : ''}`;
            item.innerHTML = `
                <i data-lucide="book-open"></i>
                <span>${name}</span>
                ${isBookmarked ? '<i data-lucide="bookmark" class="bm-icon"></i>' : ''}
            `;
            item.onclick = () => this.loadChapter(name, data[name]);
            list.appendChild(item);
        });
        lucide.createIcons();
    },

    renderDashboard(data) {
        const view = document.getElementById('chapterContent');
        if (!view) return;
        view.style.opacity = '0';
        this.activeChapter = null;
        localStorage.removeItem('lastChapter');

        setTimeout(() => {
            let cardsHtml = '';
            this.chapters.forEach(name => {
                const info = data[name];
                const isBookmarked = this.bookmarks.includes(`${this.currentClass}-${this.currentSubject}-${name}`);
                cardsHtml += `
                    <div class="chapter-card animate-slide-up"
                        onclick="app.loadChapter('${name}', app.cache['${this.currentClass}-${this.currentSubject}']['${name}'])">
                        <div class="card-top">
                            <div class="card-icon"><i data-lucide="book-open"></i></div>
                            <button class="bm-btn ${isBookmarked ? 'active' : ''}"
                                onclick="event.stopPropagation(); app.toggleBookmark('${name}')"
                                title="Bookmark">
                                <i data-lucide="bookmark"></i>
                            </button>
                        </div>
                        <h3>${name}</h3>
                        <p>${info.description || 'Access full study materials, formulas, and roadmap for this chapter.'}</p>
                        <div class="card-footer">
                            <span class="card-difficulty">${info.difficulty || 'Medium'}</span>
                            <span class="card-action">Read Now <i data-lucide="chevron-right"></i></span>
                        </div>
                    </div>
                `;
            });

            view.innerHTML = `
                <div class="dashboard-view">
                    <div class="badge">Subject Overview</div>
                    <h1>${this.currentSubject} — Class ${this.currentClass}</h1>
                    <p class="lead">Select a chapter below to begin your focused study session.</p>
                    <div class="dashboard-grid">${cardsHtml}</div>
                </div>
            `;
            view.style.opacity = '1';
            window.scrollTo(0, 0);
            const titleElOv = document.getElementById('activeChapterTitle');
            if (titleElOv) titleElOv.innerText = 'Overview';
            document.querySelectorAll('.chapter-item').forEach(item => item.classList.remove('active'));
            lucide.createIcons();
        }, 300);
    },

    loadChapter(name, info) {
        if (!info) { console.error('Missing info for chapter:', name); return; }
        this.activeChapter = name;
        localStorage.setItem('lastChapter', name);
        const titleEl = document.getElementById('activeChapterTitle');
        if (titleEl) titleEl.innerText = name;

        document.querySelectorAll('.chapter-item').forEach(item => {
            item.classList.toggle('active', item.querySelector('span')?.innerText.trim() === name);
        });

        const view = document.getElementById('chapterContent');
        if (!view) return;
        view.style.opacity = '0';

        const idx = this.chapters.indexOf(name);
        const prevChapter = idx > 0 ? this.chapters[idx - 1] : null;
        const nextChapter = idx < this.chapters.length - 1 ? this.chapters[idx + 1] : null;
        const isBookmarked = this.bookmarks.includes(`${this.currentClass}-${this.currentSubject}-${name}`);
        const cacheKey = `${this.currentClass}-${this.currentSubject}`;

        const makeSection = (title, icon, content) => `
            <section class="reading-section">
                <div class="section-header">
                    <h2><i data-lucide="${icon}"></i> ${title}</h2>
                </div>
                <div class="section-body">${content}</div>
            </section>`;

        const makeList = (arr) => arr && arr.length
            ? `<ul>${arr.map(i => `<li>${i}</li>`).join('')}</ul>` : '';

        const makeFormulaGrid = (arr) => arr && arr.length
            ? `<div class="formula-grid">${arr.map(f => `<div class="formula-card">${f}</div>`).join('')}</div>` : '';

        const makeDefinitionList = (arr) => arr && arr.length
            ? `<dl class="definition-list">${arr.map(d =>
                typeof d === 'object'
                    ? `<dt>${d.term}</dt><dd>${d.definition}</dd>`
                    : `<dt></dt><dd>${d}</dd>`
            ).join('')}</dl>` : '';

        let sections = '';

        if (info.introduction) sections += makeSection('Introduction', 'info', `<p>${info.introduction}</p>`);
        if (info.theory) sections += makeSection('Detailed Theory', 'book-open', `<div class="theory-text">${info.theory}</div>`);
        if (info.roadmap) sections += makeSection('Learning Roadmap', 'map', makeList(info.roadmap));
        if (info.topics) sections += makeSection('Core Topics', 'list-checks', makeList(info.topics));
        if (info.concepts) sections += makeSection('Key Concepts & Definitions', 'lightbulb', makeDefinitionList(info.concepts));
        if (info.laws) sections += makeSection('Laws & Principles', 'scale', makeList(info.laws));
        if (info.formulas) sections += makeSection('Key Formulas', 'sigma', makeFormulaGrid(info.formulas));
        if (info.derivations) sections += makeSection('Derivations', 'pencil', makeList(info.derivations));
        if (info.diagrams) sections += makeSection('Important Diagrams', 'image', makeList(info.diagrams));
        if (info.applications) sections += makeSection('Real-Life Applications', 'rocket', makeList(info.applications));
        if (info.examples) sections += makeSection('Solved Examples', 'calculator', makeList(info.examples));
        if (info.common_mistakes) sections += makeSection('Common Mistakes', 'alert-triangle', makeList(info.common_mistakes));
        if (info.exam_topics) sections += makeSection('Exam-Oriented Topics', 'target', makeList(info.exam_topics));
        if (info.faqs) sections += makeSection('Frequently Asked Questions', 'help-circle', makeList(info.faqs));
        if (info.revision) sections += makeSection('Quick Revision Points', 'zap', makeList(info.revision));
        if (info.importance) sections += makeSection('Exam Importance', 'star', `<p class="importance-note">${info.importance}</p>`);
        if (info.summary) sections += makeSection('Summary', 'check-circle', `<p>${info.summary}</p>`);

        setTimeout(() => {
            view.innerHTML = `
                <div class="animate-slide-up chapter-full">
                    <div class="chapter-meta">
                        <span class="difficulty-tag">${info.difficulty || 'Medium'} Difficulty</span>
                        <button class="bm-btn ${isBookmarked ? 'active' : ''}"
                            onclick="app.toggleBookmark('${name}')" id="bmBtn" title="Bookmark this chapter">
                            <i data-lucide="bookmark"></i>
                            <span>${isBookmarked ? 'Bookmarked' : 'Bookmark'}</span>
                        </button>
                    </div>
                    <h1>${name}</h1>
                    <p class="lead">${info.description || ''}</p>
                    ${sections}
                    <div class="chapter-nav-btns">
                        ${prevChapter
                    ? `<button class="nav-btn" onclick="app.loadChapter('${prevChapter}', app.cache['${cacheKey}']['${prevChapter}'])">
                                <i data-lucide="arrow-left"></i> ${prevChapter}
                              </button>`
                    : '<div></div>'}
                        ${nextChapter
                    ? `<button class="nav-btn nav-btn-next" onclick="app.loadChapter('${nextChapter}', app.cache['${cacheKey}']['${nextChapter}'])">
                                ${nextChapter} <i data-lucide="arrow-right"></i>
                              </button>`
                    : '<div></div>'}
                    </div>
                </div>
            `;
            view.style.opacity = '1';
            lucide.createIcons();
        }, 300);

        if (window.innerWidth <= 992) this.toggleSidebar(false);
    },

    toggleSection(header) {
        // collapse removed — sections always open
    },

    toggleBookmark(chapterName) {
        const key = `${this.currentClass}-${this.currentSubject}-${chapterName}`;
        if (this.bookmarks.includes(key)) {
            this.bookmarks = this.bookmarks.filter(b => b !== key);
        } else {
            this.bookmarks.push(key);
        }
        localStorage.setItem('bookmarks', JSON.stringify(this.bookmarks));

        // Update bookmark button in view
        const bmBtn = document.getElementById('bmBtn');
        if (bmBtn) {
            const isNowBookmarked = this.bookmarks.includes(key);
            bmBtn.classList.toggle('active', isNowBookmarked);
            bmBtn.querySelector('span').innerText = isNowBookmarked ? 'Bookmarked' : 'Bookmark';
        }

        // Re-render sidebar to update bookmark icons
        const cacheKey = `${this.currentClass}-${this.currentSubject}`;
        const data = this.cache[cacheKey];
        if (data) this.renderSidebar(data);
    },

    searchChapters(query) {
        query = query.toLowerCase();
        document.querySelectorAll('.chapter-item').forEach(item => {
            const text = item.querySelector('span')?.innerText.toLowerCase() || '';
            item.style.display = text.includes(query) ? 'flex' : 'none';
        });
    },

    updateProgressBar() {
        const winScroll = document.body.scrollTop || document.documentElement.scrollTop;
        const height = document.documentElement.scrollHeight - document.documentElement.clientHeight;
        const scrolled = height ? (winScroll / height) * 100 : 0;
        const bar = document.getElementById('progressBar');
        if (bar) bar.style.width = scrolled + '%';
    },

    toggleBackToTop() {
        const btn = document.getElementById('backToTop');
        if (btn) btn.classList.toggle('hidden', window.scrollY <= 300);
    },

    scrollToTop() {
        window.scrollTo({ top: 0, behavior: 'smooth' });
    },

    toggleSidebar(force) {
        const sb = document.getElementById('sidebar');
        this.isSidebarActive = force !== undefined ? force : !this.isSidebarActive;
        sb.classList.toggle('active', this.isSidebarActive);
    },
};

document.addEventListener('DOMContentLoaded', () => app.init());

// =============================================
// SHISHIMANU AI ASSISTANT
// =============================================
const shishimanu = {
    isOpen: false,
    history: [],   // { role: 'user'|'assistant', content: string }[]
    isLoading: false,

    toggle() {
        this.isOpen = !this.isOpen;
        const panel = document.getElementById('shishimanuPanel');
        if (this.isOpen) {
            panel.classList.remove('hidden');
            // Re-trigger animation
            panel.style.animation = 'none';
            panel.offsetHeight; // force reflow
            panel.style.animation = '';
            document.getElementById('shishimanuInput').focus();
        } else {
            panel.classList.add('hidden');
        }
    },

    async send() {
        if (this.isLoading) return;
        const input = document.getElementById('shishimanuInput');
        const message = input.value.trim();
        if (!message) return;

        // Append user bubble
        this._appendBubble('user', message);
        this.history.push({ role: 'user', content: message });
        input.value = '';

        // Show typing indicator
        this._setLoading(true);

        try {
            const res = await fetch('/api/shishimanu/chat', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ message, history: this.history.slice(0, -1) })
            });
            const data = await res.json();

            if (!res.ok || data.error) {
                const errMsg = data.error || 'Something went wrong. Please try again.';
                this._appendBubble('bot', '⚠️ ' + errMsg, true);
            } else {
                const reply = data.reply;
                this._appendBubble('bot', reply);
                this.history.push({ role: 'assistant', content: reply });
                // Keep history manageable
                if (this.history.length > 24) this.history = this.history.slice(-20);
            }
        } catch (err) {
            this._appendBubble('bot', '⚠️ Network error — make sure the server is running!', true);
        } finally {
            this._setLoading(false);
        }
    },

    _appendBubble(role, text, isError = false) {
        const container = document.getElementById('shishimanuMessages');
        const msgDiv = document.createElement('div');
        msgDiv.className = `shishimanu-msg ${role}`;

        const bubble = document.createElement('div');
        bubble.className = `shishimanu-bubble${isError ? ' error-bubble' : ''}`;

        // Convert newlines and basic markdown-like bold (**text**)
        bubble.innerHTML = text
            .replace(/&/g, '&amp;')
            .replace(/</g, '&lt;')
            .replace(/>/g, '&gt;')
            .replace(/\*\*(.*?)\*\*/g, '<strong>$1</strong>')
            .replace(/\*(.*?)\*/g, '<em>$1</em>')
            .replace(/\n/g, '<br>');

        msgDiv.appendChild(bubble);
        container.appendChild(msgDiv);
        container.scrollTop = container.scrollHeight;
    },

    _setLoading(loading) {
        this.isLoading = loading;
        const typing = document.getElementById('shishimanuTyping');
        const sendBtn = document.getElementById('shishimanuSendBtn');
        if (loading) {
            typing.classList.remove('hidden');
            sendBtn.disabled = true;
            // Scroll to show typing dots
            const container = document.getElementById('shishimanuMessages');
            container.scrollTop = container.scrollHeight;
        } else {
            typing.classList.add('hidden');
            sendBtn.disabled = false;
        }
    }
};

// =============================================
// DEVELOPER INFO MODAL
// =============================================
const devModal = {
    open() {
        const overlay = document.getElementById('devModalOverlay');
        overlay.classList.remove('hidden');
        overlay.offsetHeight; // force reflow for animation
        overlay.classList.add('visible');
        lucide.createIcons();
    },
    close() {
        const overlay = document.getElementById('devModalOverlay');
        overlay.classList.remove('visible');
        setTimeout(() => overlay.classList.add('hidden'), 350);
    },
    closeOnOverlay(e) {
        if (e.target === document.getElementById('devModalOverlay')) {
            this.close();
        }
    }
};

// Close dev modal on Escape key
document.addEventListener('keydown', (e) => {
    if (e.key === 'Escape') devModal.close();
});