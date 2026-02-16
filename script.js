const app = {
    currentSubject: 'Physics',
    currentClass: '11',
    isDetailOpen: false,
    onboardingClass: null,
    onboardingSubject: null,
    cache: {}, // Simple cache to store fetched data

    init() {
        // Initialize simple theme check
        if (localStorage.getItem('theme') === 'dark') {
            document.body.setAttribute('data-theme', 'dark');
        }

        // Check if user has completed onboarding before
        const hasOnboarded = localStorage.getItem('hasOnboarded');

        // Setup listeners immediately
        this.setupEventListeners();

        if (hasOnboarded) {
            // Recover saved state
            const savedClass = localStorage.getItem('userClass');
            const savedSubject = localStorage.getItem('userSubject');

            if (savedClass) this.currentClass = savedClass;
            if (savedSubject) this.currentSubject = savedSubject;

            this.showApp();
            this.loadSubject(this.currentSubject);

            // Sync UI
            const classSelect = document.getElementById('classSelect');
            if (classSelect) classSelect.value = this.currentClass;
        } else {
            // Ensure landing page is visible if not onboarded
            document.getElementById('landingPage').style.display = 'flex';
            document.getElementById('appContainer').style.display = 'none';
        }

        // Initialize lucide icons for landing page
        lucide.createIcons();
    },

    setupEventListeners() {
        const classSelect = document.getElementById('classSelect');
        if (classSelect) {
            classSelect.addEventListener('change', (e) => {
                this.currentClass = e.target.value;
                localStorage.setItem('userClass', this.currentClass);
                this.refreshContent();
            });
        }

        // Add home/reset functionality to logo
        document.querySelector('.logo-text').parentElement.addEventListener('click', () => {
            if (confirm('Return to Home Screen?')) {
                this.resetApp();
            }
        });
    },

    resetApp() {
        localStorage.removeItem('hasOnboarded');
        location.reload();
    },

    showOnboarding() {
        document.getElementById('landingPage').style.display = 'none';
        document.getElementById('onboardingModal').classList.remove('hidden');
        lucide.createIcons();
    },

    selectClass(event, classNum) {
        this.onboardingClass = classNum;
        this.currentClass = classNum;

        // Update button states
        document.querySelectorAll('[data-class]').forEach(btn => {
            btn.classList.remove('selected');
        });
        event.target.closest('button').classList.add('selected');

        this.checkOnboardingComplete();
        lucide.createIcons();
    },

    selectOnboardingSubject(event, subject) {
        this.onboardingSubject = subject;
        this.currentSubject = subject;

        // Update button states
        document.querySelectorAll('[data-subject]').forEach(btn => {
            btn.classList.remove('selected');
        });
        event.target.closest('button').classList.add('selected');

        this.checkOnboardingComplete();
        lucide.createIcons();
    },

    checkOnboardingComplete() {
        const continueBtn = document.getElementById('continueBtn');
        if (this.onboardingClass && this.onboardingSubject) {
            continueBtn.disabled = false;
        } else {
            continueBtn.disabled = true;
        }
    },

    completeOnboarding() {
        // Hide onboarding modal
        document.getElementById('onboardingModal').classList.add('hidden');

        // Save onboarding status and preferences
        localStorage.setItem('hasOnboarded', 'true');
        localStorage.setItem('userClass', this.currentClass);
        localStorage.setItem('userSubject', this.currentSubject);

        // Show the app
        this.showApp();

        // Load selected subject
        this.loadSubject(this.currentSubject);
        this.setupEventListeners();

        // Update class select
        document.getElementById('classSelect').value = this.currentClass;

        lucide.createIcons();
    },

    showApp() {
        document.getElementById('landingPage').style.display = 'none';
        document.getElementById('appContainer').style.display = 'flex';
    },

    toggleTheme() {
        const isDark = document.body.getAttribute('data-theme') === 'dark';
        document.body.setAttribute('data-theme', isDark ? 'light' : 'dark');
        localStorage.setItem('theme', isDark ? 'light' : 'dark');
    },

    loadSubject(subject) {
        this.currentSubject = subject;

        // Update Sidebar Active State
        document.querySelectorAll('.nav-item').forEach(btn => {
            btn.classList.remove('active');
            if (btn.dataset.subject === subject) btn.classList.add('active');
        });

        // Update Header
        const titles = {
            'Physics': 'Physics Dashboard',
            'Chemistry': 'Chemistry Laboratory',
            'Maths': 'Mathematics Center'
        };
        const pageTitle = document.getElementById('pageTitle');
        if (pageTitle) {
            pageTitle.innerText = titles[subject] || 'Dashboard';
        }

        this.refreshContent();
    },

    async fetchData(classId, subject) {
        const cacheKey = `${classId}-${subject}`;
        if (this.cache[cacheKey]) {
            return this.cache[cacheKey];
        }

        try {
            const response = await fetch(`/api/class/${classId}/subject/${subject}`);
            const data = await response.json();
            this.cache[cacheKey] = data;
            return data;
        } catch (error) {
            console.error('Error fetching data:', error);
            return null;
        }
    },

    async refreshContent() {
        const grid = document.getElementById("chapterGrid");
        if (!grid) return;

        grid.innerHTML = "<div style='color: var(--gold); text-align: center; width: 100%; grid-column: 1/-1;'>Loading modules...</div>";

        const subjectData = await this.fetchData(this.currentClass, this.currentSubject);

        grid.innerHTML = "";

        if (!subjectData || Object.keys(subjectData).length === 0) {
            grid.innerHTML = "<p>No data available for this selection.</p>";
            return;
        }

        let delay = 0;
        Object.entries(subjectData).forEach(([chapterName, info]) => {
            const card = document.createElement("div");
            card.className = "card animate-fade-in";
            card.style.animationDelay = `${delay}ms`;

            // Get icon based on subject
            let iconName = 'book-open';
            if (this.currentSubject === 'Physics') iconName = 'atom';
            if (this.currentSubject === 'Chemistry') iconName = 'flask-conical';
            if (this.currentSubject === 'Maths') iconName = 'calculator';

            card.innerHTML = `
                <div class="card-icon">
                    <i data-lucide="${iconName}"></i>
                </div>
                <div>
                    <h3>${chapterName}</h3>
                    <p style="color: var(--text-secondary); font-size: 0.9rem;">
                        ${info.topics.slice(0, 2).join(", ")}...
                    </p>
                </div>
                <div class="card-meta">
                    <span class="difficulty-badge">${info.difficulty}</span>
                    <button class="explore-btn">
                        Explore Chapter <i data-lucide="arrow-right" style="width: 14px;"></i>
                    </button>
                </div>
            `;

            // Make whole card clickable
            card.addEventListener('click', () => {
                this.showDetail(chapterName, info);
            });

            // Make button clickable with stopPropagation
            const btn = card.querySelector('.explore-btn');
            if (btn) {
                btn.addEventListener('click', (e) => {
                    e.stopPropagation();
                    this.showDetail(chapterName, info);
                });
            }

            grid.appendChild(card);
            delay += 50;
        });

        lucide.createIcons();
    },

    searchChapters(query) {
        query = query.toLowerCase();
        document.querySelectorAll(".card").forEach(card => {
            const title = card.querySelector("h3").innerText.toLowerCase();
            if (title.includes(query)) {
                card.style.display = "flex";
            } else {
                card.style.display = "none";
            }
        });
    },

    showDetail(chapter, info) {
        const panel = document.getElementById("detailPanel");
        const content = document.getElementById("detailContent");

        // Generate Roadmap HTML safely
        const roadmapHTML = Array.isArray(info.roadmap)
            ? info.roadmap.map(step => `<div class="roadmap-step"><span>${step}</span></div>`).join('')
            : `<p>${info.roadmap}</p>`;

        // Add description section if it exists
        const descriptionHTML = info.description
            ? `
            <div class="section-card" style="border-left: 4px solid var(--gold);">
                <h4 class="section-title"><i data-lucide="info"></i> Chapter Overview</h4>
                <p style="color: var(--text-primary); font-size: 1.05rem; line-height: 1.6;">${info.description}</p>
            </div>`
            : '';

        content.innerHTML = `
            <div class="detail-header">
                <span style="color: var(--accent); font-weight: 600; font-size: 0.9rem; text-transform: uppercase; letter-spacing: 1px;">Chapter Detail</span>
                <h2 class="detail-title">${chapter}</h2>
                <span class="difficulty-badge" style="font-size: 0.9rem;">${info.difficulty} Difficulty</span>
            </div>

            ${descriptionHTML}

            <div class="section-card">
                <h4 class="section-title"><i data-lucide="map"></i> Learning Roadmap</h4>
                <div class="roadmap-container">
                    ${roadmapHTML}
                </div>
            </div>

            <div class="section-card">
                <h4 class="section-title"><i data-lucide="list"></i> Topics Covered</h4>
                <ul style="list-style-position: inside; color: var(--text-secondary);">
                    ${info.topics.map(t => `<li style="margin-bottom: 0.5rem;">${t}</li>`).join("")}
                </ul>
            </div>

            <div class="section-card">
                <h4 class="section-title"><i data-lucide="function-square"></i> Key Formulas</h4>
                <div class="formulas-grid">
                    ${info.formulas.map(f => `<div class="formula-box">${f}</div>`).join("")}
                </div>
            </div>

            <div class="section-card">
                <h4 class="section-title"><i data-lucide="alert-circle"></i> Exam Importance</h4>
                <p style="color: var(--text-secondary);">${info.importance}</p>
            </div>
        `;

        panel.classList.remove('hidden'); // Ensure it's visible before sliding
        // Small delay to ensure display:block applies before transition starts
        requestAnimationFrame(() => {
            panel.classList.add('active'); // Slide in
        });

        this.isDetailOpen = true;
        lucide.createIcons();
    },

    closeDetail() {
        const panel = document.getElementById("detailPanel");
        panel.classList.remove('active');
        this.isDetailOpen = false;

        // Wait for transition to finish (0.5s) before hiding
        setTimeout(() => {
            if (!this.isDetailOpen) { // Check in case it was reopened
                panel.classList.add('hidden');
            }
        }, 500);
    }
};

// Start the app
document.addEventListener('DOMContentLoaded', () => {
    app.init();
});
