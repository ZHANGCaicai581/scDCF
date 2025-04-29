"""
Check and install dependencies for scDCF package.
"""

import importlib
import subprocess
import sys
import logging

logger = logging.getLogger(__name__)

REQUIRED_PACKAGES = {
    'numpy': 'numpy',
    'pandas': 'pandas',
    'scanpy': 'scanpy',
    'anndata': 'anndata',
    'scipy': 'scipy',
    'statsmodels': 'statsmodels',
    'matplotlib': 'matplotlib',
    'seaborn': 'seaborn',
    'tqdm': 'tqdm'
}

def check_and_install_dependencies():
    """
    Check if all required dependencies are installed and install missing ones.
    """
    missing_packages = []
    
    for module_name, package_name in REQUIRED_PACKAGES.items():
        try:
            importlib.import_module(module_name)
            logger.debug(f"Package {package_name} is already installed.")
        except ImportError:
            missing_packages.append(package_name)
            logger.warning(f"Package {package_name} is not installed.")
    
    if missing_packages:
        logger.info(f"Installing missing packages: {', '.join(missing_packages)}")
        try:
            subprocess.check_call([sys.executable, "-m", "pip", "install"] + missing_packages)
            logger.info("All missing packages have been installed successfully.")
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to install missing packages: {e}")
            logger.error("Please install the missing packages manually.")
            raise ImportError("Required packages are missing and could not be installed automatically.")
    else:
        logger.info("All required packages are already installed.") 