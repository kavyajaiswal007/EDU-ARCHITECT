from flask import Flask, jsonify, send_from_directory, request
import os
try:
    import anthropic
    _anthropic_available = True
except ImportError:
    _anthropic_available = False


app = Flask(__name__)

# Data storage
jls_extract_var = "Thermodynamics"
data = {
    "11": {
        "Physics": {
            jls_extract_var: {
                "difficulty": "Hard",
                "description": "Thermodynamics is the branch of physics that deals with the relationship between heat, work, and energy. It explains how thermal energy is converted to other forms of energy and how it affects matter. The laws of thermodynamics govern the principles of energy conservation, the direction of heat flow, and the concept of entropy in physical and chemical processes.",
                "roadmap": ["Heat vs Temperature concepts", "Zeroth & First Law", "Work & Internal Energy", "Second Law & Carnot Cycle", "Entropy & Irreversibility"],
                "topics": ["Zeroth Law - Thermal Equilibrium", "First Law - Energy Conservation", "Second Law - Heat Flow Direction", "Carnot Engine - Ideal Efficiency", "Entropy - Disorder Measure"],
                "formulas": ["Q = mcΔT (Heat Transfer)", "ΔU = Q − W (First Law)", "η = 1 − T₂/T₁ (Carnot Efficiency)", "W = PΔV (Work Done)"],
                "importance": "Critical for JEE Main & Advanced. 10-12% weightage in boards. Used in engines, refrigerators."
            },
            "Kinematics": {
                "difficulty": "Medium",
                "description": "Kinematics is the study of motion of objects without considering the forces that cause the motion. It describes how objects move using concepts like displacement, velocity, and acceleration. From straight-line motion to projectiles, kinematics provides the mathematical framework to predict and analyze the trajectory and speed of moving bodies.",
                "roadmap": ["Distance vs Displacement", "Velocity & Acceleration", "Equations of Motion", "Projectile Motion", "Relative Velocity"],
                "topics": ["Uniform & Non-uniform Motion", "Average & Instantaneous Velocity", "Acceleration in 1D & 2D", "Projectile Trajectories", "Relative Motion Problems"],
                "formulas": ["v = u + at", "s = ut + ½at²", "v² = u² + 2as", "R = u²sin2θ/g (Range)", "H = u²sin²θ/2g (Max Height)"],
                "importance": "Foundation of mechanics. 15% board weightage. Essential for JEE problem-solving."
            },
            "Laws of Motion": {
                "difficulty": "Hard",
                "description": "Laws of Motion, formulated by Sir Isaac Newton, form the foundation of classical mechanics. They describe the relationship between a body and the forces acting upon it, and the body's motion in response to those forces. These three laws explain inertia, force-acceleration relationship, and action-reaction pairs that govern every physical interaction.",
                "roadmap": ["Newton's Three Laws", "Free Body Diagrams", "Friction Types", "Circular Motion", "Connected Bodies"],
                "topics": ["Inertia & Mass", "Linear Momentum", "Impulse-Momentum Theorem", "Static & Kinetic Friction", "Uniform Circular Motion"],
                "formulas": ["F = ma", "p = mv", "J = FΔt = Δp", "f = μN", "F = mv²/r (Centripetal)"],
                "importance": "Core concept with 20% board weightage. Appears in 8-10 JEE questions annually."
            },
            "Work, Energy & Power": {
                "difficulty": "Medium",
                "description": "Work, Energy & Power explores how forces transfer energy to objects. Work is done when a force moves an object, energy is the capacity to do work, and power measures the rate of doing work. This chapter connects mechanical concepts through the Work-Energy theorem and the principle of conservation of energy — one of nature's most fundamental laws.",
                "roadmap": ["Work by Constant Force", "Kinetic & Potential Energy", "Work-Energy Theorem", "Conservation of Energy", "Power & Efficiency"],
                "topics": ["Work Done Calculation", "Kinetic Energy Derivation", "Gravitational & Elastic PE", "Conservative Forces", "Average & Instantaneous Power"],
                "formulas": ["W = F·s·cosθ", "KE = ½mv²", "PE = mgh", "Power = W/t", "η = Output/Input"],
                "importance": "12-15% board weightage. Crucial for mechanics problems in competitive exams."
            },
            "Gravitation": {
                "difficulty": "Hard",
                "description": "Gravitation deals with the universal attractive force between all masses. Newton's Law of Universal Gravitation explains how every object attracts every other object with a force proportional to their masses. This chapter covers planetary motion, satellite orbits, escape velocity, and Kepler's laws that govern celestial mechanics.",
                "roadmap": ["Universal Law", "Gravitational Field", "Potential & Potential Energy", "Orbital Motion", "Kepler's Laws"],
                "topics": ["Newton's Law of Gravitation", "Gravitational Field Intensity", "Escape Velocity", "Satellite Motion", "Planetary Motion"],
                "formulas": ["F = Gm₁m₂/r²", "g = GM/R²", "Vₑ = √(2GM/R)", "T² ∝ r³ (Kepler's 3rd Law)"],
                "importance": "10% board weightage. Important for space science applications in JEE."
            },
            "Rotational Motion": {
                "difficulty": "Hard",
                "description": "Rotational Motion studies objects spinning around a fixed axis. It extends Newton's laws to rotating bodies using torque, moment of inertia, and angular momentum. This chapter explains why a spinning top stays upright, how figure skaters spin faster by pulling their arms in, and the physics behind wheels, gears, and rolling objects.",
                "roadmap": ["Angular Kinematics", "Torque & Moment of Inertia", "Angular Momentum", "Rolling Motion", "Equilibrium of Rigid Bodies"],
                "topics": ["Angular Displacement & Velocity", "Torque Calculation", "Parallel & Perpendicular Axis Theorem", "Conservation of Angular Momentum", "Pure Rolling Condition"],
                "formulas": ["τ = r × F", "L = Iω", "I = Σmr²", "KE = ½Iω² + ½mv²", "α = τ/I"],
                "importance": "15% board weightage. Complex problems in JEE Advanced with high scoring potential."
            },
            "Simple Harmonic Motion": {
                "difficulty": "Medium",
                "description": "Simple Harmonic Motion (SHM) is a type of periodic oscillatory motion where the restoring force is directly proportional to the displacement. It describes the motion of pendulums, springs, and vibrating strings. SHM is foundational for understanding waves, sound, and alternating current in advanced physics.",
                "roadmap": ["Periodic Motion Basics", "SHM Characteristics", "Spring-Mass System", "Simple Pendulum", "Energy in SHM"],
                "topics": ["Displacement, Velocity, Acceleration in SHM", "Time Period & Frequency", "Spring Constant", "Pendulum Motion", "Total Energy Conservation"],
                "formulas": ["x = A sin(ωt + φ)", "T = 2π√(m/k)", "T = 2π√(l/g)", "E = ½kA²"],
                "importance": "10-12% board weightage. Appears in JEE with wave motion integration."
            },
            "Waves": {
                "difficulty": "Medium",
                "description": "Waves chapter explores how energy propagates through a medium without the physical transfer of matter. It covers transverse and longitudinal waves, wave equations, superposition principle, interference, beats, and standing waves. Understanding waves is essential for sound, light, and electromagnetic radiation studies.",
                "roadmap": ["Wave Motion Types", "Wave Equation", "Speed of Waves", "Superposition Principle", "Standing Waves"],
                "topics": ["Transverse & Longitudinal Waves", "Wave Speed, Frequency, Wavelength", "Interference & Beats", "Stationary Waves", "Resonance"],
                "formulas": ["v = fλ", "y = A sin(kx - ωt)", "f_beat = |f₁ - f₂|", "v = √(T/μ) (String)"],
                "importance": "8-10% board weightage. Foundation for sound and optics in Class 12."
            }
        },
        "Chemistry": {
            "Chemical Bonding": {
                "difficulty": "Hard",
                "description": "Chemical Bonding explains how atoms combine to form molecules and compounds by sharing or transferring electrons. It covers ionic, covalent, and metallic bonds, along with theories like VSEPR and Molecular Orbital Theory that predict molecular shapes and properties.",
                "roadmap": ["Octet Rule & Exceptions", "Ionic vs Covalent Bonding", "Lewis Structures", "VSEPR Theory", "Hybridization & MO Theory"],
                "topics": ["Electronegativity & Bond Polarity", "Lewis Dot Structures", "VSEPR - Molecular Geometry", "sp, sp², sp³ Hybridization", "Molecular Orbital Theory"],
                "formulas": ["Bond order = (Nb - Na)/2", "% Ionic Character = 16|Δχ| + 3.5(Δχ)²"],
                "importance": "Foundation for organic & inorganic. 12-15% board weightage. Critical for JEE."
            },
            "Atomic Structure": {
                "difficulty": "Medium",
                "description": "Atomic Structure explores the internal architecture of atoms — from subatomic particles (protons, neutrons, electrons) to quantum mechanical models. It explains electron energy levels, orbitals, and the rules governing electron configuration that determine an element's chemical behavior.",
                "roadmap": ["Thomson & Rutherford Models", "Bohr's Atomic Model", "Quantum Mechanical Model", "Quantum Numbers", "Electronic Configuration"],
                "topics": ["Subatomic Particles", "Atomic Spectra", "Energy Levels", "Orbitals s,p,d,f", "Aufbau, Pauli, Hund's Rules"],
                "formulas": ["E = hν", "λ = h/mv (de Broglie)", "1/λ = R(1/n₁² - 1/n₂²)", "ΔE = E₂ - E₁"],
                "importance": "Core concept. 10-12% board weightage. Basics for entire chemistry."
            },
            "Periodic Properties": {
                "difficulty": "Medium",
                "description": "Periodic Properties examines how element characteristics change across the periodic table. It covers trends in atomic radius, ionization energy, electron affinity, and electronegativity — the predictable patterns that allow chemists to anticipate how elements will behave in chemical reactions.",
                "roadmap": ["Periodic Table Organization", "Atomic Radius Trends", "Ionization & Electron Affinity", "Electronegativity", "Metallic Character"],
                "topics": ["Groups & Periods", "Effective Nuclear Charge", "Shielding Effect", "Ionization Energy Trends", "Electronegativity Scales"],
                "formulas": ["IE₁ < IE₂ < IE₃...", "Zeff = Z - σ"],
                "importance": "8-10% board weightage. Essential for predicting chemical behavior."
            },
            "Thermodynamics": {
                "difficulty": "Hard",
                "description": "Chemical Thermodynamics studies energy changes during chemical reactions — whether reactions release or absorb heat. It introduces enthalpy, entropy, and Gibbs free energy to predict whether a reaction will occur spontaneously and how much energy it involves.",
                "roadmap": ["System & Surroundings", "Internal Energy & Enthalpy", "First Law", "Hess's Law", "Entropy & Gibbs Energy"],
                "topics": ["Heat & Work", "Enthalpy of Formation", "Bond Enthalpy", "Spontaneity", "Gibb's Free Energy"],
                "formulas": ["ΔU = q + w", "ΔH = ΔU + PΔV", "ΔG = ΔH - TΔS", "ΔHrxn = ΣΔHproducts - ΣΔHreactants"],
                "importance": "12-15% board weightage. High scoring in JEE Advanced."
            },
            "States of Matter": {
                "difficulty": "Easy",
                "description": "States of Matter explains the behavior of gases, liquids, and solids under different conditions of temperature and pressure. It covers gas laws, the kinetic molecular theory, intermolecular forces, and properties like surface tension and viscosity.",
                "roadmap": ["Gas Laws", "Kinetic Theory", "Liquefaction", "Liquid State Properties"],
                "topics": ["Boyle's, Charles's Law", "Ideal Gas Equation", "Graham's Law", "Van der Waals Equation", "Surface Tension, Viscosity"],
                "formulas": ["PV = nRT", "r₁/r₂ = √(M₂/M₁)", "(P + an²/V²)(V - nb) = nRT"],
                "importance": "10% board weightage. Numerical-based questions in JEE."
            },
            "Equilibrium": {
                "difficulty": "Hard",
                "description": "Equilibrium explores the state where forward and reverse reactions occur at equal rates. It covers Le Chatelier's Principle, equilibrium constants Kc and Kp, ionic equilibrium including acids, bases, pH calculations, and buffer solutions.",
                "roadmap": ["Chemical Equilibrium Concept", "Equilibrium Constant", "Le Chatelier's Principle", "Ionic Equilibrium", "Buffer Solutions"],
                "topics": ["Kc & Kp Relations", "Reaction Quotient Q", "Acid-Base Equilibrium", "pH & pOH", "Common Ion Effect"],
                "formulas": ["Kp = Kc(RT)^Δn", "pH = -log[H⁺]", "Henderson-Hasselbalch: pH = pKa + log([A⁻]/[HA])"],
                "importance": "15-18% board weightage. Maximum JEE questions from this chapter."
            },
            "Redox Reactions": {
                "difficulty": "Medium",
                "description": "Redox Reactions involve the transfer of electrons between chemical species — one substance is oxidized (loses electrons) while another is reduced (gains electrons). This chapter builds the foundation for electrochemistry and explains corrosion, batteries, and metabolic reactions.",
                "roadmap": ["Oxidation Number Rules", "Balancing Redox Reactions", "Electrochemical Series", "Types of Redox"],
                "topics": ["Oxidation & Reduction", "Half-Reaction Method", "Oxidizing & Reducing Agents", "Disproportionation"],
                "formulas": ["ON change method", "Electrons lost = Electrons gained"],
                "importance": "8-10% board weightage. Connects to electrochemistry in Class 12."
            },
            "Organic Chemistry Basics": {
                "difficulty": "Medium",
                "description": "Organic Chemistry Basics introduces the chemistry of carbon compounds — the foundation of life itself. It covers IUPAC naming, functional groups, isomerism, reaction mechanisms (SN1, SN2, electrophilic, nucleophilic), and electronic effects like inductive and mesomeric effects.",
                "roadmap": ["IUPAC Nomenclature", "Isomerism Types", "Reaction Mechanisms", "GOC - General Organic Chemistry"],
                "topics": ["Functional Groups", "Structural & Stereoisomerism", "Inductive, Mesomeric Effects", "Hyperconjugation", "Electrophile & Nucleophile"],
                "formulas": ["+I, -I effects", "+M, -M effects", "Markovnikov's Rule"],
                "importance": "Foundation for Class 12 organic. 12% board weightage."
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
            "Electrostatics": {
                "difficulty": "Hard",
                "roadmap": ["Electric Charge Properties", "Coulomb's Law", "Electric Field & Field Lines", "Gauss's Law", "Electric Potential & Capacitance"],
                "topics": ["Charge Quantization & Conservation", "Superposition Principle", "Electric Field Intensity", "Electric Flux", "Capacitors in Series & Parallel"],
                "formulas": ["F = kq₁q₂/r² = (1/4πε₀)q₁q₂/r²", "E = F/q", "V = W/q", "C = Q/V", "C = ε₀A/d"],
                "importance": "High weightage 15-18% in boards. Core JEE topic with numerical problems."
            },
            "Current Electricity": {
                "difficulty": "Medium",
                "roadmap": ["Electric Current & Drift Velocity", "Ohm's Law", "Resistance & Resistivity", "Kirchhoff's Laws", "Wheatstone Bridge"],
                "topics": ["Current Density", "Temperature Dependence", "Series & Parallel Combinations", "Meter Bridge", "Potentiometer"],
                "formulas": ["I = Q/t", "V = IR", "R = ρl/A", "Rs = R₁ + R₂...", "1/Rp = 1/R₁ + 1/R₂..."],
                "importance": "12-15% board weightage. Circuit problems frequent in JEE."
            },
            "Magnetic Effects of Current": {
                "difficulty": "Hard",
                "roadmap": ["Biot-Savart Law", "Ampere's Law", "Force on Current Carrying Conductor", "Torque on Coil", "Moving Coil Galvanometer"],
                "topics": ["Magnetic Field due to Straight Wire", "Circular Loop", "Solenoid", "Lorentz Force", "Magnetic Dipole Moment"],
                "formulas": ["B = μ₀I/2πr", "B = μ₀NI/l (Solenoid)", "F = BIlsinθ", "τ = NIABsinθ"],
                "importance": "15% board weightage. Conceptual JEE questions with high difficulty."
            },
            "Electromagnetic Induction": {
                "difficulty": "Hard",
                "roadmap": ["Faraday's Law", "Lenz's Law", "Motional EMF", "Self & Mutual Inductance", "AC Generators"],
                "topics": ["Magnetic Flux", "Induced EMF & Current", "Eddy Currents", "Energy in Inductor", "Transformers"],
                "formulas": ["ε = -dΦ/dt", "ε = Blv", "L = Φ/I", "M = Φ₂/I₁", "ε = NABω sinωt"],
                "importance": "18-20% board weightage. Maximum JEE Main & Advanced problems."
            },
            "Optics": {
                "difficulty": "Medium",
                "roadmap": ["Reflection Laws", "Refraction & Snell's Law", "Lenses & Mirrors", "Wave Optics", "Optical Instruments"],
                "topics": ["Mirror Formula", "Lens Formula", "Total Internal Reflection", "Interference & Diffraction", "Young's Double Slit"],
                "formulas": ["1/f = 1/v + 1/u", "μ = sin i / sin r", "μ = c/v", "β = λD/d (Fringe Width)"],
                "importance": "15% board weightage. Ray & wave optics both important for JEE."
            },
            "Modern Physics": {
                "difficulty": "Hard",
                "roadmap": ["Photoelectric Effect", "Matter Waves", "Bohr's Model", "Nuclear Physics", "Radioactivity"],
                "topics": ["Work Function", "Einstein's Equation", "de Broglie Wavelength", "Mass-Energy Equivalence", "Half-Life & Decay"],
                "formulas": ["E = hν = hc/λ", "KE = hν - φ", "λ = h/p", "E = mc²", "N = N₀(½)^(t/T)"],
                "importance": "12-15% board weightage. Conceptual clarity needed for JEE Advanced."
            },
            "Semiconductor Electronics": {
                "difficulty": "Medium",
                "roadmap": ["Intrinsic & Extrinsic Semiconductors", "PN Junction Diode", "Diode Applications", "Transistors", "Logic Gates"],
                "topics": ["Doping - n-type & p-type", "Forward & Reverse Bias", "Rectifiers", "Transistor as Amplifier", "Digital Electronics"],
                "formulas": ["I = I₀(e^(eV/kT) - 1)", "β = Ic/Ib", "α = Ic/Ie"],
                "importance": "10% board weightage. Application-based JEE questions."
            },
            "Communication Systems": {
                "difficulty": "Easy",
                "roadmap": ["Communication Elements", "Bandwidth & Range", "Modulation", "AM & FM", "Internet & Mobile"],
                "topics": ["Transmitter & Receiver", "Signal Propagation", "Amplitude Modulation", "Frequency Modulation", "Satellite Communication"],
                "formulas": ["Modulation Index m = Am/Ac", "Bandwidth = 2fm"],
                "importance": "5-8% board weightage. Theory-based easy marks."
            }
        },
        "Chemistry": {
            "Electrochemistry": {
                "difficulty": "Medium",
                "roadmap": ["Redox Reactions Review", "Galvanic & Electrolytic Cells", "Nernst Equation", "Conductance", "Batteries & Fuel Cells"],
                "topics": ["Cell Notation", "Standard Electrode Potential", "EMF Calculation", "Kohlrausch's Law", "Corrosion Prevention"],
                "formulas": ["E°cell = E°cathode - E°anode", "E = E° - (0.0591/n) logQ", "ΔG = -nFE", "Λm = κ/c"],
                "importance": "15% board weightage. Numerical problems frequent in JEE Main."
            },
            "Chemical Kinetics": {
                "difficulty": "Hard",
                "roadmap": ["Rate of Reaction", "Rate Law & Order", "Integrated Rate Equations", "Arrhenius Equation", "Collision Theory"],
                "topics": ["Zero, First, Second Order Reactions", "Half-life Period", "Temperature Dependence", "Activation Energy", "Catalysis"],
                "formulas": ["Rate = k[A]^m[B]^n", "t1/2 = 0.693/k (First Order)", "k = Ae^(-Ea/RT)"],
                "importance": "12-15% board weightage. Graph-based JEE questions."
            },
            "d & f Block Elements": {
                "difficulty": "Medium",
                "roadmap": ["Electronic Configuration", "Oxidation States", "Color & Magnetism", "Lanthanoids", "Actinoids"],
                "topics": ["Transition Elements Properties", "K2Cr2O7 & KMnO4", "Coordination Compounds", "Inner Transition Elements"],
                "formulas": ["Magnetic Moment μ = √[n(n+2)] BM"],
                "importance": "10% board weightage. Inorganic chemistry JEE favorite."
            },
            "Coordination Compounds": {
                "difficulty": "Hard",
                "roadmap": ["Werner's Theory", "IUPAC Nomenclature", "Isomerism", "Bonding Theories", "Crystal Field Theory"],
                "topics": ["Ligands & Coordination Number", "Geometrical & Optical Isomerism", "VBT", "CFT - Splitting", "Stability Constants"],
                "formulas": ["EAN = Z - oxidation number + 2×CN", "CFSE Calculation"],
                "importance": "12% board weightage. Complex naming & structure in JEE."
            },
            "Aldehydes, Ketones & Carboxylic Acids": {
                "difficulty": "Hard",
                "roadmap": ["Nomenclature", "Carbonyl Group Reactions", "Aldol Condensation", "Cannizzaro Reaction", "Carboxylic Acid Properties"],
                "topics": ["Nucleophilic Addition", "Oxidation & Reduction", "Acidity of Carboxylic Acids", "Esterification", "Derivatives"],
                "formulas": ["Tollens' Test", "Fehling's Test", "Iodoform Test"],
                "importance": "18-20% board weightage. Maximum organic chemistry marks."
            },
            "Amines": {
                "difficulty": "Medium",
                "roadmap": ["Classification", "Preparation Methods", "Basicity", "Diazonium Salts", "Aromatic Amines"],
                "topics": ["Gabriel Phthalimide Synthesis", "Hoffmann Bromamide", "Hinsberg's Test", "Coupling Reactions", "Aniline Reactions"],
                "formulas": ["Basicity Order: Aliphatic > Aromatic", "pKb trends"],
                "importance": "10-12% board weightage. Name reactions important for JEE."
            },
            "Polymers": {
                "difficulty": "Easy",
                "roadmap": ["Classification", "Addition Polymers", "Condensation Polymers", "Biodegradable Polymers", "Applications"],
                "topics": ["Homopolymers & Copolymers", "Nylon, Terylene, Bakelite", "Natural & Synthetic Rubber", "Plastics"],
                "formulas": ["Degree of Polymerization n = Molecular weight of polymer / Molecular weight of monomer"],
                "importance": "6-8% board weightage. Easy theory-based marks."
            },
            "Biomolecules": {
                "difficulty": "Medium",
                "roadmap": ["Carbohydrates", "Proteins", "Nucleic Acids", "Vitamins", "Enzymes"],
                "topics": ["Monosaccharides, Oligosaccharides", "Amino Acids & Peptides", "DNA & RNA Structure", "Vitamin Classification"],
                "formulas": ["Glucose Structure", "Peptide Bond Formation"],
                "importance": "8-10% board weightage. Structural understanding for JEE."
            },
            "Chemistry in Everyday Life": {
                "difficulty": "Easy",
                "roadmap": ["Drugs & Medicines", "Chemotherapy", "Antibiotics", "Detergents", "Food Preservatives"],
                "topics": ["Analgesics, Antipyretics", "Antiseptics vs Disinfectants", "Tranquilizers", "Antacids", "Artificial Sweeteners"],
                "formulas": ["No major formulas - Conceptual"],
                "importance": "5% board weightage. Easy scoring one-liners for JEE."
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
# SHISHIMANU AI ASSISTANT — Claude-Powered
# =========================================

@app.route('/api/shishimanu/chat', methods=['POST'])
def shishimanu_chat():
    if not _anthropic_available:
        return jsonify({'error': 'anthropic package not installed. Run: pip install anthropic'}), 500

    api_key = os.environ.get('ANTHROPIC_API_KEY')
    if not api_key:
        return jsonify({'error': 'ANTHROPIC_API_KEY not set. Please set it as an environment variable.'}), 500

    body = request.get_json(force=True, silent=True) or {}
    user_message = body.get('message', '').strip()
    history = body.get('history', [])  # list of {role, content} dicts

    if not user_message:
        return jsonify({'error': 'No message provided'}), 400

    # Build messages array
    messages = []
    for h in history[-10:]:  # keep last 10 turns for context
        role = h.get('role', 'user')
        content = h.get('content', '')
        if role in ('user', 'assistant') and content:
            messages.append({'role': role, 'content': content})
    messages.append({'role': 'user', 'content': user_message})

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
        client = anthropic.Anthropic(api_key=api_key)
        response = client.messages.create(
            model='claude-3-5-haiku-20241022',
            max_tokens=1024,
            system=system_prompt,
            messages=messages
        )
        reply = response.content[0].text
        return jsonify({'reply': reply})
    except anthropic.AuthenticationError:
        return jsonify({'error': 'Invalid API key. Please check your ANTHROPIC_API_KEY.'}), 401
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/favicon.ico')
def favicon():
    return '', 204

@app.route('/<path:path>')
def serve_static(path):
    return send_from_directory('.', path)

if __name__ == '__main__':
    app.run(debug=True, port=3000)
