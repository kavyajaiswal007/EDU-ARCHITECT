const CHAPTER_DATA = {
    "11": {
        "Physics": {
            "Thermodynamics": {
                "difficulty": "Hard",
                "description": "Thermodynamics is the branch of physics that deals with the relationship between heat, work, and energy. It explains how thermal energy is converted to other forms of energy and how it affects matter.",
                "roadmap": ["Heat vs Temperature concepts", "Zeroth & First Law", "Work & Internal Energy", "Second Law & Carnot Cycle", "Entropy & Irreversibility"],
                "topics": ["Zeroth Law - Thermal Equilibrium", "First Law - Energy Conservation", "Second Law - Heat Flow Direction", "Carnot Engine - Ideal Efficiency", "Entropy - Disorder Measure"],
                "formulas": ["Q = mcΔT (Heat Transfer)", "ΔU = Q − W (First Law)", "η = 1 − T₂/T₁ (Carnot Efficiency)", "W = PΔV (Work Done)"],
                "importance": "Critical for JEE Main & Advanced. 10-12% weightage in boards. Used in engines, refrigerators."
            },
            "Kinematics": {
                "difficulty": "Medium",
                "description": "Kinematics is the study of motion of objects without considering the forces that cause the motion. It describes how objects move using concepts like displacement, velocity, and acceleration.",
                "roadmap": ["Distance vs Displacement", "Velocity & Acceleration", "Equations of Motion", "Projectile Motion", "Relative Velocity"],
                "topics": ["Uniform & Non-uniform Motion", "Average & Instantaneous Velocity", "Acceleration in 1D & 2D", "Projectile Trajectories", "Relative Motion Problems"],
                "formulas": ["v = u + at", "s = ut + ½at²", "v² = u² + 2as", "R = u²sin2θ/g (Range)", "H = u²sin²θ/2g (Max Height)"],
                "importance": "Foundation of mechanics. 15% board weightage. Essential for JEE problem-solving."
            },
            "Laws of Motion": {
                "difficulty": "Hard",
                "description": "Laws of Motion, formulated by Sir Isaac Newton, form the foundation of classical mechanics. They describe the relationship between a body and the forces acting upon it.",
                "roadmap": ["Newton's Three Laws", "Free Body Diagrams", "Friction Types", "Circular Motion", "Connected Bodies"],
                "topics": ["Inertia & Mass", "Linear Momentum", "Impulse-Momentum Theorem", "Static & Kinetic Friction", "Uniform Circular Motion"],
                "formulas": ["F = ma", "p = mv", "J = FΔt = Δp", "f = μN", "F = mv²/r (Centripetal)"],
                "importance": "Core concept with 20% board weightage. Appears in 8-10 JEE questions annually."
            },
            "Work, Energy & Power": {
                "difficulty": "Medium",
                "description": "Work, Energy & Power explores how forces transfer energy to objects. Work is done when a force moves an object, energy is the capacity to do work, and power measures the rate of doing work.",
                "roadmap": ["Work by Constant Force", "Kinetic & Potential Energy", "Work-Energy Theorem", "Conservation of Energy", "Power & Efficiency"],
                "topics": ["Work Done Calculation", "Kinetic Energy Derivation", "Gravitational & Elastic PE", "Conservative Forces", "Average & Instantaneous Power"],
                "formulas": ["W = F·s·cosθ", "KE = ½mv²", "PE = mgh", "Power = W/t", "η = Output/Input"],
                "importance": "12-15% board weightage. Crucial for mechanics problems in competitive exams."
            },
            "Gravitation": {
                "difficulty": "Hard",
                "description": "Gravitation deals with the universal attractive force between all masses. Newton's Law of Universal Gravitation explains how every object attracts every other with a force proportional to their masses.",
                "roadmap": ["Universal Law", "Gravitational Field", "Potential & Potential Energy", "Orbital Motion", "Kepler's Laws"],
                "topics": ["Newton's Law of Gravitation", "Gravitational Field Intensity", "Escape Velocity", "Satellite Motion", "Planetary Motion"],
                "formulas": ["F = Gm₁m₂/r²", "g = GM/R²", "Vₑ = √(2GM/R)", "T² ∝ r³ (Kepler's 3rd Law)"],
                "importance": "10% board weightage. Important for space science applications in JEE."
            },
            "Rotational Motion": {
                "difficulty": "Hard",
                "description": "Rotational Motion studies objects spinning around a fixed axis. It extends Newton's laws to rotating bodies using torque, moment of inertia, and angular momentum.",
                "roadmap": ["Angular Kinematics", "Torque & Moment of Inertia", "Angular Momentum", "Rolling Motion", "Equilibrium of Rigid Bodies"],
                "topics": ["Angular Displacement & Velocity", "Torque Calculation", "Parallel & Perpendicular Axis Theorem", "Conservation of Angular Momentum", "Pure Rolling Condition"],
                "formulas": ["τ = r × F", "L = Iω", "I = Σmr²", "KE = ½Iω² + ½mv²", "α = τ/I"],
                "importance": "15% board weightage. Complex problems in JEE Advanced with high scoring potential."
            },
            "Simple Harmonic Motion": {
                "difficulty": "Medium",
                "description": "Simple Harmonic Motion (SHM) is a type of periodic oscillatory motion where the restoring force is directly proportional to the displacement.",
                "roadmap": ["Periodic Motion Basics", "SHM Characteristics", "Spring-Mass System", "Simple Pendulum", "Energy in SHM"],
                "topics": ["Displacement, Velocity, Acceleration in SHM", "Time Period & Frequency", "Spring Constant", "Pendulum Motion", "Total Energy Conservation"],
                "formulas": ["x = A sin(ωt + φ)", "T = 2π√(m/k)", "T = 2π√(l/g)", "E = ½kA²"],
                "importance": "10-12% board weightage. Appears in JEE with wave motion integration."
            },
            "Waves": {
                "difficulty": "Medium",
                "description": "Waves chapter explores how energy propagates through a medium without the physical transfer of matter. It covers transverse and longitudinal waves, wave equations, and standing waves.",
                "roadmap": ["Wave Motion Types", "Wave Equation", "Speed of Waves", "Superposition Principle", "Standing Waves"],
                "topics": ["Transverse & Longitudinal Waves", "Wave Speed, Frequency, Wavelength", "Interference & Beats", "Stationary Waves", "Resonance"],
                "formulas": ["v = fλ", "y = A sin(kx - ωt)", "f_beat = |f₁ - f₂|", "v = √(T/μ) (String)"],
                "importance": "8-10% board weightage. Foundation for sound and optics in Class 12."
            }
        },
        "Chemistry": {
            "Chemical Bonding": {
                "difficulty": "Hard",
                "description": "Chemical Bonding explains how atoms combine to form molecules and compounds by sharing or transferring electrons. It covers ionic, covalent, and metallic bonds, along with VSEPR and Molecular Orbital Theory.",
                "roadmap": ["Octet Rule & Exceptions", "Ionic vs Covalent Bonding", "Lewis Structures", "VSEPR Theory", "Hybridization & MO Theory"],
                "topics": ["Electronegativity & Bond Polarity", "Lewis Dot Structures", "VSEPR - Molecular Geometry", "sp, sp², sp³ Hybridization", "Molecular Orbital Theory"],
                "formulas": ["Bond order = (Nb - Na)/2", "% Ionic Character = 16|Δχ| + 3.5(Δχ)²"],
                "importance": "Foundation for organic & inorganic. 12-15% board weightage. Critical for JEE."
            },
            "Atomic Structure": {
                "difficulty": "Medium",
                "description": "Atomic Structure explores the internal architecture of atoms — from subatomic particles to quantum mechanical models. It explains electron energy levels, orbitals, and electron configuration.",
                "roadmap": ["Thomson & Rutherford Models", "Bohr's Atomic Model", "Quantum Mechanical Model", "Quantum Numbers", "Electronic Configuration"],
                "topics": ["Subatomic Particles", "Atomic Spectra", "Energy Levels", "Orbitals s,p,d,f", "Aufbau, Pauli, Hund's Rules"],
                "formulas": ["E = hν", "λ = h/mv (de Broglie)", "1/λ = R(1/n₁² - 1/n₂²)", "ΔE = E₂ - E₁"],
                "importance": "Core concept. 10-12% board weightage. Basics for entire chemistry."
            },
            "Periodic Properties": {
                "difficulty": "Medium",
                "description": "Periodic Properties examines how element characteristics change across the periodic table — atomic radius, ionization energy, electron affinity, and electronegativity trends.",
                "roadmap": ["Periodic Table Organization", "Atomic Radius Trends", "Ionization & Electron Affinity", "Electronegativity", "Metallic Character"],
                "topics": ["Groups & Periods", "Effective Nuclear Charge", "Shielding Effect", "Ionization Energy Trends", "Electronegativity Scales"],
                "formulas": ["IE₁ < IE₂ < IE₃...", "Zeff = Z - σ"],
                "importance": "8-10% board weightage. Essential for predicting chemical behavior."
            },
            "Thermodynamics": {
                "difficulty": "Hard",
                "description": "Chemical Thermodynamics studies energy changes during chemical reactions — enthalpy, entropy, and Gibbs free energy to predict whether a reaction will occur spontaneously.",
                "roadmap": ["System & Surroundings", "Internal Energy & Enthalpy", "First Law", "Hess's Law", "Entropy & Gibbs Energy"],
                "topics": ["Heat & Work", "Enthalpy of Formation", "Bond Enthalpy", "Spontaneity", "Gibb's Free Energy"],
                "formulas": ["ΔU = q + w", "ΔH = ΔU + PΔV", "ΔG = ΔH - TΔS", "ΔHrxn = ΣΔHproducts - ΣΔHreactants"],
                "importance": "12-15% board weightage. High scoring in JEE Advanced."
            },
            "States of Matter": {
                "difficulty": "Easy",
                "description": "States of Matter explains the behavior of gases, liquids, and solids under different conditions. It covers gas laws, kinetic molecular theory, and intermolecular forces.",
                "roadmap": ["Gas Laws", "Kinetic Theory", "Liquefaction", "Liquid State Properties"],
                "topics": ["Boyle's, Charles's Law", "Ideal Gas Equation", "Graham's Law", "Van der Waals Equation", "Surface Tension, Viscosity"],
                "formulas": ["PV = nRT", "r₁/r₂ = √(M₂/M₁)", "(P + an²/V²)(V - nb) = nRT"],
                "importance": "10% board weightage. Numerical-based questions in JEE."
            },
            "Equilibrium": {
                "difficulty": "Hard",
                "description": "Equilibrium explores the state where forward and reverse reactions occur at equal rates. Covers Le Chatelier's Principle, equilibrium constants, ionic equilibrium, pH, and buffer solutions.",
                "roadmap": ["Chemical Equilibrium Concept", "Equilibrium Constant", "Le Chatelier's Principle", "Ionic Equilibrium", "Buffer Solutions"],
                "topics": ["Kc & Kp Relations", "Reaction Quotient Q", "Acid-Base Equilibrium", "pH & pOH", "Common Ion Effect"],
                "formulas": ["Kp = Kc(RT)^Δn", "pH = -log[H⁺]", "Henderson-Hasselbalch: pH = pKa + log([A⁻]/[HA])"],
                "importance": "15-18% board weightage. Maximum JEE questions from this chapter."
            },
            "Redox Reactions": {
                "difficulty": "Medium",
                "description": "Redox Reactions involve the transfer of electrons between chemical species. Builds the foundation for electrochemistry and explains corrosion, batteries, and metabolic reactions.",
                "roadmap": ["Oxidation Number Rules", "Balancing Redox Reactions", "Electrochemical Series", "Types of Redox"],
                "topics": ["Oxidation & Reduction", "Half-Reaction Method", "Oxidizing & Reducing Agents", "Disproportionation"],
                "formulas": ["ON change method", "Electrons lost = Electrons gained"],
                "importance": "8-10% board weightage. Connects to electrochemistry in Class 12."
            },
            "Organic Chemistry Basics": {
                "difficulty": "Medium",
                "description": "Organic Chemistry Basics introduces the chemistry of carbon compounds. Covers IUPAC naming, functional groups, isomerism, and reaction mechanisms like inductive and mesomeric effects.",
                "roadmap": ["IUPAC Nomenclature", "Isomerism Types", "Reaction Mechanisms", "GOC - General Organic Chemistry"],
                "topics": ["Functional Groups", "Structural & Stereoisomerism", "Inductive, Mesomeric Effects", "Hyperconjugation", "Electrophile & Nucleophile"],
                "formulas": ["+I, -I effects", "+M, -M effects", "Markovnikov's Rule"],
                "importance": "Foundation for Class 12 organic. 12% board weightage."
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
            "Electrostatics": {
                "difficulty": "Hard",
                "description": "Electrostatics covers electric charge properties, Coulomb's Law, electric field, Gauss's Law, electric potential, and capacitance — the foundation of all electrical engineering.",
                "roadmap": ["Electric Charge Properties", "Coulomb's Law", "Electric Field & Field Lines", "Gauss's Law", "Electric Potential & Capacitance"],
                "topics": ["Charge Quantization & Conservation", "Superposition Principle", "Electric Field Intensity", "Electric Flux", "Capacitors in Series & Parallel"],
                "formulas": ["F = kq₁q₂/r²", "E = F/q", "V = W/q", "C = Q/V", "C = ε₀A/d"],
                "importance": "High weightage 15-18% in boards. Core JEE topic with numerical problems."
            },
            "Current Electricity": {
                "difficulty": "Medium",
                "description": "Current Electricity covers electric current, Ohm's Law, resistance & resistivity, Kirchhoff's Laws, Wheatstone Bridge, and the Potentiometer.",
                "roadmap": ["Electric Current & Drift Velocity", "Ohm's Law", "Resistance & Resistivity", "Kirchhoff's Laws", "Wheatstone Bridge"],
                "topics": ["Current Density", "Temperature Dependence", "Series & Parallel Combinations", "Meter Bridge", "Potentiometer"],
                "formulas": ["I = Q/t", "V = IR", "R = ρl/A", "Rs = R₁ + R₂...", "1/Rp = 1/R₁ + 1/R₂..."],
                "importance": "12-15% board weightage. Circuit problems frequent in JEE."
            },
            "Magnetic Effects of Current": {
                "difficulty": "Hard",
                "description": "Magnetic Effects of Current covers Biot-Savart Law, Ampere's Law, force on current-carrying conductors, torque on coils, and the Moving Coil Galvanometer.",
                "roadmap": ["Biot-Savart Law", "Ampere's Law", "Force on Current Carrying Conductor", "Torque on Coil", "Moving Coil Galvanometer"],
                "topics": ["Magnetic Field due to Straight Wire", "Circular Loop", "Solenoid", "Lorentz Force", "Magnetic Dipole Moment"],
                "formulas": ["B = μ₀I/2πr", "B = μ₀NI/l (Solenoid)", "F = BIlsinθ", "τ = NIABsinθ"],
                "importance": "15% board weightage. Conceptual JEE questions with high difficulty."
            },
            "Electromagnetic Induction": {
                "difficulty": "Hard",
                "description": "Electromagnetic Induction covers Faraday's Law, Lenz's Law, Motional EMF, Self & Mutual Inductance, and AC Generators.",
                "roadmap": ["Faraday's Law", "Lenz's Law", "Motional EMF", "Self & Mutual Inductance", "AC Generators"],
                "topics": ["Magnetic Flux", "Induced EMF & Current", "Eddy Currents", "Energy in Inductor", "Transformers"],
                "formulas": ["ε = -dΦ/dt", "ε = Blv", "L = Φ/I", "M = Φ₂/I₁", "ε = NABω sinωt"],
                "importance": "18-20% board weightage. Maximum JEE Main & Advanced problems."
            },
            "Optics": {
                "difficulty": "Medium",
                "description": "Optics covers reflection, refraction, lenses & mirrors, wave optics, and optical instruments like microscopes and telescopes.",
                "roadmap": ["Reflection Laws", "Refraction & Snell's Law", "Lenses & Mirrors", "Wave Optics", "Optical Instruments"],
                "topics": ["Mirror Formula", "Lens Formula", "Total Internal Reflection", "Interference & Diffraction", "Young's Double Slit"],
                "formulas": ["1/f = 1/v + 1/u", "μ = sin i / sin r", "μ = c/v", "β = λD/d (Fringe Width)"],
                "importance": "15% board weightage. Ray & wave optics both important for JEE."
            },
            "Modern Physics": {
                "difficulty": "Hard",
                "description": "Modern Physics covers the Photoelectric Effect, Matter Waves, Bohr's Model, Nuclear Physics, and Radioactivity.",
                "roadmap": ["Photoelectric Effect", "Matter Waves", "Bohr's Model", "Nuclear Physics", "Radioactivity"],
                "topics": ["Work Function", "Einstein's Equation", "de Broglie Wavelength", "Mass-Energy Equivalence", "Half-Life & Decay"],
                "formulas": ["E = hν = hc/λ", "KE = hν - φ", "λ = h/p", "E = mc²", "N = N₀(½)^(t/T)"],
                "importance": "12-15% board weightage. Conceptual clarity needed for JEE Advanced."
            },
            "Semiconductor Electronics": {
                "difficulty": "Medium",
                "description": "Semiconductor Electronics covers intrinsic & extrinsic semiconductors, PN junction diode, rectifiers, transistors, and logic gates.",
                "roadmap": ["Intrinsic & Extrinsic Semiconductors", "PN Junction Diode", "Diode Applications", "Transistors", "Logic Gates"],
                "topics": ["Doping - n-type & p-type", "Forward & Reverse Bias", "Rectifiers", "Transistor as Amplifier", "Digital Electronics"],
                "formulas": ["I = I₀(e^(eV/kT) - 1)", "β = Ic/Ib", "α = Ic/Ie"],
                "importance": "10% board weightage. Application-based JEE questions."
            },
            "Communication Systems": {
                "difficulty": "Easy",
                "description": "Communication Systems covers elements of communication, bandwidth, modulation, AM & FM, and internet & mobile communication.",
                "roadmap": ["Communication Elements", "Bandwidth & Range", "Modulation", "AM & FM", "Internet & Mobile"],
                "topics": ["Transmitter & Receiver", "Signal Propagation", "Amplitude Modulation", "Frequency Modulation", "Satellite Communication"],
                "formulas": ["Modulation Index m = Am/Ac", "Bandwidth = 2fm"],
                "importance": "5-8% board weightage. Theory-based easy marks."
            }
        },
        "Chemistry": {
            "Electrochemistry": {
                "difficulty": "Medium",
                "description": "Electrochemistry covers galvanic & electrolytic cells, the Nernst Equation, conductance, batteries, and fuel cells.",
                "roadmap": ["Redox Reactions Review", "Galvanic & Electrolytic Cells", "Nernst Equation", "Conductance", "Batteries & Fuel Cells"],
                "topics": ["Cell Notation", "Standard Electrode Potential", "EMF Calculation", "Kohlrausch's Law", "Corrosion Prevention"],
                "formulas": ["E°cell = E°cathode - E°anode", "E = E° - (0.0591/n) logQ", "ΔG = -nFE", "Λm = κ/c"],
                "importance": "15% board weightage. Numerical problems frequent in JEE Main."
            },
            "Chemical Kinetics": {
                "difficulty": "Hard",
                "description": "Chemical Kinetics studies rate of reaction, rate law, integrated rate equations, Arrhenius Equation, and collision theory.",
                "roadmap": ["Rate of Reaction", "Rate Law & Order", "Integrated Rate Equations", "Arrhenius Equation", "Collision Theory"],
                "topics": ["Zero, First, Second Order Reactions", "Half-life Period", "Temperature Dependence", "Activation Energy", "Catalysis"],
                "formulas": ["Rate = k[A]^m[B]^n", "t1/2 = 0.693/k (First Order)", "k = Ae^(-Ea/RT)"],
                "importance": "12-15% board weightage. Graph-based JEE questions."
            },
            "d & f Block Elements": {
                "difficulty": "Medium",
                "description": "d & f Block Elements covers electronic configuration, oxidation states, color & magnetism of transition elements, lanthanoids and actinoids.",
                "roadmap": ["Electronic Configuration", "Oxidation States", "Color & Magnetism", "Lanthanoids", "Actinoids"],
                "topics": ["Transition Elements Properties", "K2Cr2O7 & KMnO4", "Coordination Compounds", "Inner Transition Elements"],
                "formulas": ["Magnetic Moment μ = √[n(n+2)] BM"],
                "importance": "10% board weightage. Inorganic chemistry JEE favorite."
            },
            "Coordination Compounds": {
                "difficulty": "Hard",
                "description": "Coordination Compounds covers Werner's Theory, IUPAC nomenclature, isomerism, VBT, and Crystal Field Theory.",
                "roadmap": ["Werner's Theory", "IUPAC Nomenclature", "Isomerism", "Bonding Theories", "Crystal Field Theory"],
                "topics": ["Ligands & Coordination Number", "Geometrical & Optical Isomerism", "VBT", "CFT - Splitting", "Stability Constants"],
                "formulas": ["EAN = Z - oxidation number + 2×CN", "CFSE Calculation"],
                "importance": "12% board weightage. Complex naming & structure in JEE."
            },
            "Aldehydes, Ketones & Carboxylic Acids": {
                "difficulty": "Hard",
                "description": "Covers carbonyl group reactions — aldol condensation, Cannizzaro reaction, nucleophilic addition, oxidation & reduction, and carboxylic acid derivatives.",
                "roadmap": ["Nomenclature", "Carbonyl Group Reactions", "Aldol Condensation", "Cannizzaro Reaction", "Carboxylic Acid Properties"],
                "topics": ["Nucleophilic Addition", "Oxidation & Reduction", "Acidity of Carboxylic Acids", "Esterification", "Derivatives"],
                "formulas": ["Tollens' Test", "Fehling's Test", "Iodoform Test"],
                "importance": "18-20% board weightage. Maximum organic chemistry marks."
            },
            "Amines": {
                "difficulty": "Medium",
                "description": "Amines covers classification, preparation methods, basicity, diazonium salts, and aromatic amines including aniline reactions.",
                "roadmap": ["Classification", "Preparation Methods", "Basicity", "Diazonium Salts", "Aromatic Amines"],
                "topics": ["Gabriel Phthalimide Synthesis", "Hoffmann Bromamide", "Hinsberg's Test", "Coupling Reactions", "Aniline Reactions"],
                "formulas": ["Basicity Order: Aliphatic > Aromatic", "pKb trends"],
                "importance": "10-12% board weightage. Name reactions important for JEE."
            },
            "Polymers": {
                "difficulty": "Easy",
                "description": "Polymers covers classification, addition & condensation polymers, biodegradable polymers, and their applications.",
                "roadmap": ["Classification", "Addition Polymers", "Condensation Polymers", "Biodegradable Polymers", "Applications"],
                "topics": ["Homopolymers & Copolymers", "Nylon, Terylene, Bakelite", "Natural & Synthetic Rubber", "Plastics"],
                "formulas": ["Degree of Polymerization n = Molecular weight of polymer / Molecular weight of monomer"],
                "importance": "6-8% board weightage. Easy theory-based marks."
            },
            "Biomolecules": {
                "difficulty": "Medium",
                "description": "Biomolecules covers carbohydrates, proteins, nucleic acids, vitamins, and enzymes — their structure and biological functions.",
                "roadmap": ["Carbohydrates", "Proteins", "Nucleic Acids", "Vitamins", "Enzymes"],
                "topics": ["Monosaccharides, Oligosaccharides", "Amino Acids & Peptides", "DNA & RNA Structure", "Vitamin Classification"],
                "formulas": ["Glucose Structure", "Peptide Bond Formation"],
                "importance": "8-10% board weightage. Structural understanding for JEE."
            },
            "Chemistry in Everyday Life": {
                "difficulty": "Easy",
                "description": "Chemistry in Everyday Life covers drugs, medicines, antibiotics, detergents, and food preservatives — their types and uses.",
                "roadmap": ["Drugs & Medicines", "Chemotherapy", "Antibiotics", "Detergents", "Food Preservatives"],
                "topics": ["Analgesics, Antipyretics", "Antiseptics vs Disinfectants", "Tranquilizers", "Antacids", "Artificial Sweeteners"],
                "formulas": ["No major formulas - Conceptual"],
                "importance": "5% board weightage. Easy scoring one-liners for JEE."
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
    },

    showApp() {
        document.getElementById('landingPage').style.display = 'none';
        document.getElementById('appContainer').style.display = 'flex';
        document.body.style.overflow = 'auto';
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
        document.getElementById('headerSubject').innerText = subject;
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
            document.getElementById('activeChapterTitle').innerText = 'Overview';
            document.querySelectorAll('.chapter-item').forEach(item => item.classList.remove('active'));
            lucide.createIcons();
        }, 300);
    },

    loadChapter(name, info) {
        if (!info) { console.error('Missing info for chapter:', name); return; }
        this.activeChapter = name;
        localStorage.setItem('lastChapter', name);
        document.getElementById('activeChapterTitle').innerText = name;

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
            <section class="reading-section collapsible">
                <div class="section-header" onclick="app.toggleSection(this)">
                    <h2><i data-lucide="${icon}"></i> ${title}</h2>
                    <i data-lucide="chevron-down" class="collapse-icon"></i>
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
        const section = header.closest('.collapsible');
        section.classList.toggle('collapsed');
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

