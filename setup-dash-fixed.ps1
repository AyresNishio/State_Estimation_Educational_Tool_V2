# State Estimation Educational Tool - Docker Startup Script (PowerShell)
# Optimized for Dash/Plotly application

Write-Host "🎓 State Estimation Educational Tool - Docker Setup" -ForegroundColor Green
Write-Host "==================================================" -ForegroundColor Green
Write-Host "🔬 Dash/Plotly Application with Multi-language Support" -ForegroundColor Cyan
Write-Host ""

# Check if Docker is installed
try {
    docker --version | Out-Null
    Write-Host "✅ Docker is installed" -ForegroundColor Green
} catch {
    Write-Host "❌ Docker is not installed. Please install Docker Desktop first." -ForegroundColor Red
    Write-Host "Visit: https://docs.docker.com/get-docker/" -ForegroundColor Yellow
    exit 1
}

# Check if Docker is running
try {
    docker info | Out-Null
    Write-Host "✅ Docker is running" -ForegroundColor Green
} catch {
    Write-Host "❌ Docker daemon is not running. Please start Docker Desktop." -ForegroundColor Red
    exit 1
}

# Function to display menu
function Show-Menu {
    Write-Host ""
    Write-Host "Choose setup option:" -ForegroundColor Cyan
    Write-Host "1) Quick Start (Recommended) - Use existing dockerfile" -ForegroundColor White
    Write-Host "2) Optimized Setup - Use optimized Dockerfile" -ForegroundColor White
    Write-Host "3) Development Mode - Docker Compose with hot reload" -ForegroundColor White
    Write-Host "4) Production Mode - Docker Compose optimized" -ForegroundColor White
    Write-Host "5) Clean up containers" -ForegroundColor White
    Write-Host "6) View application logs" -ForegroundColor White
    Write-Host "7) Exit" -ForegroundColor White
    Write-Host ""
}

# Function for quick start
function Start-QuickSetup {
    Write-Host "🚀 Quick Start - Using existing configuration..." -ForegroundColor Yellow

    if (Test-Path "dockerfile") {
        Write-Host "✅ Found existing dockerfile" -ForegroundColor Green
        docker build -t state-estimation-tool .
        docker run -d -p 8050:8050 --name state-estimation-app state-estimation-tool
        Write-Host "✅ Application is running at http://localhost:8050" -ForegroundColor Green
        Write-Host "🌍 Multilingual support: PT-BR and EN available" -ForegroundColor Cyan
    } else {
        Write-Host "❌ dockerfile not found. Please ensure you're in the project directory." -ForegroundColor Red
        exit 1
    }
}

# Function for optimized setup
function Start-OptimizedSetup {
    Write-Host "⚡ Optimized Setup - Using enhanced Dockerfile..." -ForegroundColor Yellow

    if (Test-Path "Dockerfile.optimized") {
        Copy-Item "Dockerfile.optimized" "Dockerfile"
        docker build -t state-estimation-tool-optimized .
        docker run -d -p 8050:8050 --name state-estimation-app state-estimation-tool-optimized
        Write-Host "✅ Optimized application is running at http://localhost:8050" -ForegroundColor Green
        Write-Host "🔋 Enhanced with health checks and security features" -ForegroundColor Cyan
    } else {
        Write-Host "❌ Dockerfile.optimized not found." -ForegroundColor Red
        exit 1
    }
}

# Function for development mode
function Start-DevSetup {
    Write-Host "🛠️ Development Mode - Docker Compose with hot reload..." -ForegroundColor Yellow

    if (Test-Path "docker-compose.yml") {
        docker-compose up --build -d
        Write-Host "✅ Development environment is running at http://localhost:8050" -ForegroundColor Green
        Write-Host "🔄 Hot reload enabled - changes will reflect automatically" -ForegroundColor Cyan
        Write-Host "📊 Redis cache available at localhost:6379" -ForegroundColor Cyan
    } else {
        Write-Host "❌ docker-compose.yml not found." -ForegroundColor Red
        exit 1
    }
}

# Function for production mode
function Start-ProdSetup {
    Write-Host "🏭 Production Mode - Optimized Docker Compose..." -ForegroundColor Yellow

    if ((Test-Path "docker-compose.yml") -and (Test-Path "Dockerfile.optimized")) {
        Copy-Item "Dockerfile.optimized" "Dockerfile"
        docker-compose up --build -d
        Write-Host "✅ Production environment is running at http://localhost:8050" -ForegroundColor Green
        Write-Host "⚡ Multiple workers for better performance" -ForegroundColor Cyan
        Write-Host "🔒 Security optimizations enabled" -ForegroundColor Cyan
    } else {
        Write-Host "❌ Required files not found." -ForegroundColor Red
        exit 1
    }
}

# Function to clean up
function Start-Cleanup {
    Write-Host "🧹 Cleaning up containers and images..." -ForegroundColor Yellow

    try { docker stop state-estimation-app } catch { }
    try { docker rm state-estimation-app } catch { }
    try { docker-compose down } catch { }

    try { docker rmi state-estimation-tool } catch { }
    try { docker rmi state-estimation-tool-optimized } catch { }

    Write-Host "✅ Cleanup complete" -ForegroundColor Green
}

# Function to view logs
function Show-Logs {
    Write-Host "📋 Viewing application logs..." -ForegroundColor Yellow

    $runningContainer = docker ps --filter name=state-estimation-app --quiet
    if ($runningContainer) {
        Write-Host "Container logs:" -ForegroundColor Cyan
        docker logs --tail=50 -f state-estimation-app
    } else {
        Write-Host "❌ No running containers found" -ForegroundColor Red
    }
}

# Main menu loop
do {
    Show-Menu
    $choice = Read-Host "Enter choice [1-7]"

    switch ($choice) {
        "1" {
            Start-QuickSetup
            break
        }
        "2" {
            Start-OptimizedSetup
            break
        }
        "3" {
            Start-DevSetup
            break
        }
        "4" {
            Start-ProdSetup
            break
        }
        "5" {
            Start-Cleanup
        }
        "6" {
            Show-Logs
        }
        "7" {
            Write-Host "👋 Goodbye!" -ForegroundColor Green
            exit 0
        }
        default {
            Write-Host "❌ Invalid option. Please choose 1-7." -ForegroundColor Red
        }
    }
} while ($choice -ne "7")

if ([int]$choice -lt 5) {
    Write-Host ""
    Write-Host "🎉 Setup complete!" -ForegroundColor Green
    Write-Host ""
    Write-Host "📱 Access your State Estimation Educational Tool:" -ForegroundColor Cyan
    Write-Host "   🌐 Web Interface: http://localhost:8050" -ForegroundColor White
    Write-Host "   🇵🇹 Portuguese: http://localhost:8050/dash/pt/overview" -ForegroundColor White
    Write-Host "   🇺🇸 English: http://localhost:8050/dash/en/overview" -ForegroundColor White
    Write-Host ""
    Write-Host "📋 Useful commands:" -ForegroundColor Cyan
    Write-Host "   View logs: docker logs state-estimation-app" -ForegroundColor White
    Write-Host "   Stop app: docker stop state-estimation-app" -ForegroundColor White
    Write-Host "   Restart: docker restart state-estimation-app" -ForegroundColor White
    Write-Host "   Remove: docker rm state-estimation-app" -ForegroundColor White
}
