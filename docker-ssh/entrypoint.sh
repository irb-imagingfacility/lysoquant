#!/bin/bash
# Container entrypoint: an SFTP-only sshd for the Fiji clients, and
# lysoquant-serve working on the folder they upload to, as the same
# unprivileged user.
set -euo pipefail

NAME=lysoquant
ID_UID=${LYSOQUANT_UID:-1000}
ID_GID=${LYSOQUANT_GID:-1000}
JOBS=${LYSOQUANT_JOBS:-/jobs}
PORT=${LYSOQUANT_SSH_PORT:-2222}
HOSTKEYS=/etc/lysoquant/hostkeys

# Public keys allowed to connect: the file mounted at /authorized_keys, and
# those in $LYSOQUANT_AUTHORIZED_KEYS (several separated by ;)
KEYS_FILE=/authorized_keys
if [ -d "$KEYS_FILE" ]; then
    echo "ERROR: $KEYS_FILE is a directory: the file to mount (LYSOQUANT_AUTHORIZED_KEYS_FILE) does not exist on the host" >&2
    exit 1
fi
printf '%s\n' "${LYSOQUANT_AUTHORIZED_KEYS:-}" | tr ';' '\n' \
    | sed 's/^[[:space:]]*//;s/[[:space:]]*$//' | grep -E '^(ssh-|ecdsa-|sk-)' > /etc/lysoquant/authorized_keys || true
if [ ! -s /etc/lysoquant/authorized_keys ] && ! grep -qE '^(ssh-|ecdsa-|sk-)' "$KEYS_FILE" 2> /dev/null; then
    echo "ERROR: no public key in $KEYS_FILE or LYSOQUANT_AUTHORIZED_KEYS, nobody could connect (see README.md)" >&2
    exit 1
fi
chmod 644 /etc/lysoquant/authorized_keys

getent group "$ID_GID" > /dev/null || groupadd -g "$ID_GID" $NAME
id -u $NAME > /dev/null 2>&1 || useradd -u "$ID_UID" -g "$ID_GID" -m -s /bin/bash -p '*' $NAME

mkdir -p "$HOSTKEYS" "$JOBS"
for type in ed25519 ecdsa rsa; do
    [ -f "$HOSTKEYS/ssh_host_${type}_key" ] || ssh-keygen -q -t "$type" -N '' -f "$HOSTKEYS/ssh_host_${type}_key"
done
chown "$ID_UID:$ID_GID" "$JOBS"

/usr/sbin/sshd -D -e -f /etc/lysoquant/sshd_config -o "Port=$PORT" &
sshd=$!
setpriv --reuid="$ID_UID" --regid="$ID_GID" --init-groups env HOME=/home/$NAME lysoquant-serve "$@" &
serve=$!

# docker stop: let the running jobs finish, then stop sshd
stop() {
    kill -TERM "$serve" 2> /dev/null || true
    wait "$serve" || true
    kill -TERM "$sshd" 2> /dev/null || true
    exit 0
}
trap stop TERM INT

# If either process dies, stop the container so that Docker restarts it
status=0
wait -n || status=$?
kill -TERM "$serve" "$sshd" 2> /dev/null || true
wait || true
exit "$((status == 0 ? 1 : status))"
