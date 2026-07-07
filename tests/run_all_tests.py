#!/usr/bin/env python
import os
import logging
import subprocess
import sys

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def run_tests():
    """Run maintained package verification tests."""
    tests = [
        ([sys.executable, os.path.join("tests", "test_pipeline_smoke.py")], "End-to-end pipeline smoke test"),
        ([sys.executable, os.path.join("tests", "test_organization.py")], "Organized output regression test"),
    ]
    
    results = []
    env = os.environ.copy()
    env.setdefault("MPLCONFIGDIR", "/private/tmp/mpl")
    env.setdefault("XDG_CACHE_HOME", "/private/tmp/cache")
    env.setdefault("NUMBA_CACHE_DIR", "/private/tmp/numba")
    env["PYTHONPATH"] = REPO_ROOT + os.pathsep + env.get("PYTHONPATH", "")

    for cmd, desc in tests:
        logging.info(f"Running: {desc}")
        print(f"\n{'='*50}")
        print(f"TEST: {desc}")
        print(f"{'='*50}")
        
        result = subprocess.run(cmd, env=env)
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
