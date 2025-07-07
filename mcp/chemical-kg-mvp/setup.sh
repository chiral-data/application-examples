#!/bin/bash

# Chemical Knowledge Graph MVP Setup Script

set -e  # Exit on any error

echo "Chemical Knowledge Graph MVP - Setup Script"
echo "==========================================="

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Check if Docker is installed
check_docker() {
    print_status "Checking Docker installation..."
    if ! command -v docker &> /dev/null; then
        print_error "Docker is not installed. Please install Docker first."
        echo "Visit: https://docs.docker.com/get-docker/"
        exit 1
    fi
    print_success "Docker is installed"
}

# Check if Docker Compose is installed
check_docker_compose() {
    print_status "Checking Docker Compose installation..."
    if ! command -v docker-compose &> /dev/null && ! docker compose version &> /dev/null; then
        print_error "Docker Compose is not installed. Please install Docker Compose first."
        echo "Visit: https://docs.docker.com/compose/install/"
        exit 1
    fi
    print_success "Docker Compose is available"
}

# Create necessary directories
create_directories() {
    print_status "Creating project directories..."
    
    directories=("temp" "data" "scripts" "vectorstore")
    
    for dir in "${directories[@]}"; do
        if [ ! -d "$dir" ]; then
            mkdir -p "$dir"
            print_success "Created directory: $dir"
        else
            print_status "Directory already exists: $dir"
        fi
    done
}

# Create .env file if it doesn't exist
create_env_file() {
    print_status "Setting up environment configuration..."
    
    if [ ! -f ".env" ]; then
        if [ -f ".env.example" ]; then
            cp .env.example .env
            print_success "Created .env file from .env.example"
        else
            # Create basic .env file
            cat > .env << EOF
# Ollama Configuration
OLLAMA_MODEL=mistral
EMBEDDING_MODEL=nomic-embed-text

# DECIMER Configuration
DECIMER_API_URL=https://www.decimer.ai/process_image

# Application Settings
CHUNK_SIZE=1000
CHUNK_OVERLAP=200

# Development Settings
PYTHONUNBUFFERED=1
EOF
            print_success "Created basic .env file"
        fi
        print_warning "Please review and customize the .env file for your needs"
    else
        print_status ".env file already exists"
    fi
}

# Check system requirements
check_system_requirements() {
    print_status "Checking system requirements..."
    
    # Check available memory
    if command -v free &> /dev/null; then
        total_mem=$(free -g | awk '/^Mem:/{print $2}')
        if [ "$total_mem" -lt 8 ]; then
            print_warning "System has less than 8GB RAM. Ollama models may run slowly."
        else
            print_success "System has sufficient memory ($total_mem GB)"
        fi
    fi
    
    # Check disk space
    if command -v df &> /dev/null; then
        available_space=$(df -BG . | awk 'NR==2{print $4}' | sed 's/G//')
        if [ "$available_space" -lt 10 ]; then
            print_warning "Less than 10GB disk space available. Consider freeing up space."
        else
            print_success "Sufficient disk space available"
        fi
    fi
}

# Download Ollama models
setup_ollama_models() {
    print_status "Setting up Ollama models..."
    print_status "This will be handled by the Docker compose setup"
    print_status "Models will be downloaded on first startup"
}

# Validate docker-compose.yml
validate_compose() {
    print_status "Validating Docker Compose configuration..."
    
    if docker-compose config &> /dev/null || docker compose config &> /dev/null; then
        print_success "Docker Compose configuration is valid"
    else
        print_error "Docker Compose configuration has issues"
        exit 1
    fi
}

# Build Docker images
build_images() {
    print_status "Building Docker images..."
    
    if command -v docker-compose &> /dev/null; then
        docker-compose build
    else
        docker compose build
    fi
    
    print_success "Docker images built successfully"
}

# Final setup instructions
print_instructions() {
    echo ""
    echo "=========================================="
    print_success "Setup completed successfully!"
    echo "=========================================="
    echo ""
    echo "Next steps:"
    echo "1. Review the .env file and customize if needed"
    echo "2. Start the application with:"
    echo "   docker-compose up -d"
    echo ""
    echo "3. Access the application at:"
    echo "   http://localhost:8501"
    echo ""
    echo "4. Check service status with:"
    echo "   docker-compose ps"
    echo ""
    echo "5. View logs with:"
    echo "   docker-compose logs -f"
    echo ""
    echo "For more information, see README.md"
    echo ""
}

# Main execution
main() {
    print_status "Starting Chemical Knowledge Graph MVP setup..."
    
    check_docker
    check_docker_compose
    check_system_requirements
    create_directories
    create_env_file
    validate_compose
    build_images
    setup_ollama_models
    
    print_instructions
}

# Run main function
main "$@"