# EDU-ARCHITECT: Advanced Learning Dashboard

## Overview
EDU-ARCHITECT is a premium educational web application designed for Class 11 and 12 students. It provides a structured learning path for Physics, Chemistry, and Mathematics, complete with detailed roadmaps, topic breakdowns, formulas, and exam importance insights.

The project features a sleek, "Black & Gold" luxury aesthetic and a robust client-server architecture.

## Features
*   **Dynamic Data Backend**: A Python Flask server (`app.py`) manages all educational content, ensuring scalability and easy updates.
*   **Interactive Frontend**: A responsive HTML/JS interface that fetches data dynamically without page reloads.
*   **Premium Design**: Custom CSS variables, glassmorphism effects, and smooth animations found in `style.css`.
*   **Onboarding Flow**: A dedicated "Get Started" experience to personalize the dashboard for the student's class and subject.
*   **Chapter Deep-Dives**: Interactive "Explore Chapter" cards that reveal detailed roadmaps, formulas, and key topics.
*   **Search**: Real-time filtering of chapter cards.

## Project Structure
*   **`app.py`**: The Flask backend. Stores the dictionary of educational data and defines API routes to serve it.
*   **`index.html`**: The main entry point. Contains the structure for the Landing Page, Onboarding Modal, and Main Dashboard.
*   **`script.js`**: Handles logic. Fetches data from the backend, generates HTML cards, manages state (Class/Subject), and handles UI interactions.
*   **`style.css`**: Contains all styling rules, animations, and the comprehensive "Dark & Gold" theme system.
*   **`start_server.sh`**: A helper script to set up the virtual environment and launch the server.

## How to Run
The server is currently running. You can access the application at:
**[http://localhost:3000](http://localhost:3000)**

If you need to restart it manually:
1.  Open a terminal in the project folder.
2.  Run: `./start_server.sh`

## Recent Updates
*   **Backend Migration**: Moved hardcoded data from JS to Python for better architecture.
*   **UI Polish**: Upgraded the Landing Page with a premium radial gradient background and glowing animations.
*   **Bug Fixes**: Resolved an issue where the "Chapter Detail" panel would not appear due to conflicting CSS classes.
*   **Favicon**: Added a handler to prevent 404 errors for the site icon.

## Technologies
*   **Frontend**: HTML5, CSS3, JavaScript (ES6+), Lucide Icons
*   **Backend**: Python 3, Flask framework
