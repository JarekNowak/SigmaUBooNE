# Source AFTER setup_env.sh: makes GENIE use genie::HypQELTracePXSec for QEL-CC-LAMBDA.
HYPFIX_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export GXMLPATH="$HYPFIX_DIR/xml:$GXMLPATH"          # ModelConfiguration.xml + HypQELTracePXSec.xml override
export LD_LIBRARY_PATH="$HYPFIX_DIR/lib:$LD_LIBRARY_PATH"   # libGHypFix.so + rootmap (ROOT autoloads the class)
