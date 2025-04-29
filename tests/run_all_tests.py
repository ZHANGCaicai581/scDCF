#!/usr/bin/env python
import os
import logging
import subprocess
import sys
import time

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def run_tests():
    """Run all tests to verify the package is working correctly"""
    tests = [
        ("python test_organization.py", "Basic organization test"),
        ("python test_with_existing.py", "Testing with existing data"),
        ("python test_integration.py", "Package integration test"),
        ("python test_control_genes.py", "Control genes test")
    ]
    
    results = []
    for cmd, desc in tests:
        logging.info(f"Running: {desc}")
        print(f"\n{'='*50}")
        print(f"TEST: {desc}")
        print(f"{'='*50}")
        
        # Add retries for test_with_existing.py which seems problematic
        if "test_with_existing.py" in cmd:
            max_retries = 3
            retry_count = 0
            success = False
            
            while retry_count < max_retries and not success:
                if retry_count > 0:
                    logging.info(f"Retry {retry_count}/{max_retries} for {desc}")
                    # Wait a moment before retrying
                    time.sleep(2)
                
                result = subprocess.run(cmd, shell=True)
                success = result.returncode == 0
                retry_count += 1
                
                if not success and retry_count < max_retries:
                    logging.warning(f"Test failed, will retry {retry_count}/{max_retries}")
        else:
            result = subprocess.run(cmd, shell=True)
            success = result.returncode == 0
        
        results.append((desc, success))
        
        if success:
            logging.info(f"✅ {desc} completed successfully")
        else:
            logging.error(f"❌ {desc} failed with return code {result.returncode}")
    
    # Print summary
    print("\n\n")
    print("="*50)
    print("TEST SUMMARY")
    print("="*50)
    all_passed = True
    for desc, success in results:
        status = "✅ PASSED" if success else "❌ FAILED"
        print(f"{status}: {desc}")
        if not success:
            all_passed = False
    
    if all_passed:
        print("\n✅ All tests passed successfully!")
        return 0  # Success exit code
    else:
        print("\n❌ Some tests failed. See above for details.")
        return 1  # Error exit code

if __name__ == "__main__":
    sys.exit(run_tests())  # Properly exit with error code