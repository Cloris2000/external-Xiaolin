#!/usr/bin/env bash
# One-time setup: passwordless rsync from Trillium -> CAMH (xzhou@192.197.205.74).
#
# IMPORTANT: CAMH login via 192.197.205.74 is rssh-restricted
# (only scp / sftp / rsync — no interactive shell). So ssh-copy-id cannot work.
#
# Run on Trillium:
#   bash setup_camh_ssh_key.sh
#
# Then on CAMH (normal login, e.g. dev03 / Cursor), install the uploaded pubkey:
#   cat ~/trillium_camh.pub >> ~/.ssh/authorized_keys
#   chmod 600 ~/.ssh/authorized_keys
#   rm ~/trillium_camh.pub

set -euo pipefail

SRC_HOST="${SRC_HOST:-xzhou@192.197.205.74}"
KEY="${HOME}/.ssh/id_ed25519_camh"
REMOTE_PUB="${REMOTE_PUB:-trillium_camh.pub}"

mkdir -p "${HOME}/.ssh"
chmod 700 "${HOME}/.ssh"

if [[ ! -f "${KEY}" ]]; then
  echo "Generating SSH key: ${KEY}"
  ssh-keygen -t ed25519 -f "${KEY}" -N "" -C "trillium-to-camh-$(whoami)"
else
  echo "Key already exists: ${KEY}"
fi

# Optional convenience host alias
CFG="${HOME}/.ssh/config"
if ! grep -q "Host camh-xfer" "${CFG}" 2>/dev/null; then
  cat >> "${CFG}" <<EOF

Host camh-xfer
  HostName 192.197.205.74
  User xzhou
  IdentityFile ${KEY}
  IdentitiesOnly yes
  ServerAliveInterval 60
  ServerAliveCountMax 3
EOF
  chmod 600 "${CFG}"
  echo "Added SSH config alias: camh-xfer"
fi

echo
echo "Your public key:"
echo "----------------------------------------------------------------"
cat "${KEY}.pub"
echo "----------------------------------------------------------------"
echo
echo "Uploading pubkey to ${SRC_HOST}:~/${REMOTE_PUB}"
echo "(enter CAMH password once — scp is allowed under rssh)"
scp -i "${KEY}" "${KEY}.pub" "${SRC_HOST}:${REMOTE_PUB}"

echo
echo "=============================================================="
echo "Next step — on CAMH (full shell login, NOT via 192.197.205.74):"
echo "=============================================================="
echo "  cat ~/${REMOTE_PUB} >> ~/.ssh/authorized_keys"
echo "  chmod 600 ~/.ssh/authorized_keys"
echo "  rm ~/${REMOTE_PUB}"
echo
echo "Then back on Trillium, test with rsync (not ssh shell):"
echo "  rsync -avP -e \"ssh -i ${KEY} -o IdentitiesOnly=yes -o BatchMode=yes\" \\"
echo "    ${SRC_HOST}:/external/rprshnas01/external_data/rosmap/metadata/ROSmaster.rds /tmp/"
echo
echo "If that works without a password, run:  bash transfer_rosmap.sh"
