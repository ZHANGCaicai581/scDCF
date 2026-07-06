"""
Dependency checks for the scDCF package.
"""

import importlib
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
    Check whether required dependencies are importable.

    Packages should be installed through the package manager at environment setup
    time, not from inside the CLI at runtime.
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
        missing_str = ", ".join(missing_packages)
        logger.error(f"Missing required packages: {missing_str}")
        raise ImportError(
            "Missing required scDCF dependencies: "
            f"{missing_str}. Install them with your environment manager before running scDCF."
        )
    else:
        logger.info("All required packages are already installed.")
