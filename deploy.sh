#!/bin/bash

# deploy.sh - Helper script for deployment

echo "🚀 Proteomics Dashboard Deployment Helper"
echo "=========================================="

# Check if remote is already set
if git remote get-url origin >/dev/null 2>&1; then
    echo "✅ Git remote 'origin' already configured"
    git remote -v
else
    echo "❌ No git remote configured"
    echo "Please run: git remote add origin https://github.com/yourusername/proteomics-dashboard.git"
    echo "Replace 'yourusername' with your actual GitHub username"
    exit 1
fi

echo ""
echo "📊 Repository Status:"
git status --short

echo ""
echo "📁 Files to be deployed:"
echo "- deploy_app.py (main application)"
echo "- results/ (all pre-computed analysis results)"
echo "- utils/ (helper modules)"
echo "- requirements.txt (Python dependencies)"

echo ""
echo "📋 Pre-deployment checklist:"
echo "✅ Data folder excluded from git"
echo "✅ Test files removed"
echo "✅ Import issues fixed"
echo "✅ Dependencies optimized"

echo ""
echo "🔄 Pushing to GitHub..."
git push -u origin main

if [ $? -eq 0 ]; then
    echo ""
    echo "🎉 Successfully pushed to GitHub!"
    echo ""
    echo "📌 Next Steps:"
    echo "1. Go to https://share.streamlit.io"
    echo "2. Connect your GitHub account"
    echo "3. Select your repository"
    echo "4. Set main file to: deploy_app.py"
    echo "5. Click Deploy!"
    echo ""
    echo "🌐 Your app will be available at: https://yourappname.streamlit.app"
else
    echo "❌ Push failed. Please check your GitHub configuration."
fi
