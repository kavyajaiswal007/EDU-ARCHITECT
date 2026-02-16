from flask import Flask, jsonify, send_from_directory
import os

app = Flask(__name__)

# Data storage
data = {
    "11": {
        "Physics": {
            "Thermodynamics": {
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

@app.route('/')
def home():
    return send_from_directory('.', 'index.html')

@app.route('/<path:path>')
def serve_static(path):
    return send_from_directory('.', path)

@app.route('/api/data')
def get_data():
    return jsonify(data)

@app.route('/api/class/<class_id>/subject/<subject>')
def get_subject_data(class_id, subject):
    return jsonify(data.get(class_id, {}).get(subject, {}))

@app.route('/favicon.ico')
def favicon():
    return '', 204

if __name__ == '__main__':
    app.run(debug=True, port=3000)
