@echo off
cd /d D:\NGS_Project\VCF_Generation\final_dashboard
start http://localhost:5000
wsl bash -lc "cd /mnt/d/NGS_Project/VCF_Generation/final_dashboard && python3 app.py"
pause