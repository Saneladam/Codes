In Bash, `trap 'command' SIGNAL` executes `'command'` when the shell receives that specific signal (SIGINT for Ctrl+C; EXIT is a special builtin called at script termination). Python replaces the default interpreter traceback with a custom function registered through `signal.signal`.

### bash-clean.sh
```bash
cleanup() {
    echo "Cleaning up resources..."
    rm -f /tmp/lockfile  # example action
}

trap cleanup INT TERM EXIT

while true; do sleep 1; done
```
