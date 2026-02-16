
#!/bin/bash
if [ ! -d "venv" ]; then
    python3 -m venv venv
    venv/bin/pip install flask
fi

venv/bin/python app.py
