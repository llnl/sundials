from repligit import fetch_pack, ls_remote, send_pack
from os import environ as env

def main():
    src_remote_url = env.get("SRC_REMOTE_URL", None)
    dest_remote_url = env.get("DEST_REMOTE_URL", None)
    target_ref = env.get("CI_COMMIT_REF_NAME", None)

    # Authentication credentials
    # Note: Only provide credentials when authentication is required
    # src_username = "<username>"  # Uncomment if source repo requires auth
    # src_password = "<token>"     # Uncomment if source repo requires auth
    dest_username = env.get("GITLAB_USER_NAME", "")  # For destination repo write access
    dest_password = env.get("GITLAB_TOKEN", "")  # For destination repo write access

    # List references from source repository (without authentication)
    gh_refs = ls_remote(src_remote_url)

    # List references from destination repository (with authentication)
    gl_refs = ls_remote(dest_remote_url, username=dest_username, password=dest_password)

    want_sha = gh_refs[target_ref]
    have_shas = set(gl_refs.values())

    from_sha = gl_refs.get(target_ref) or ("0" * 40)

    if want_sha in have_shas:
        print("Everything is up to date")
        return

    # Fetch the packfile from source repository
    packfile = fetch_pack(
        src_remote_url,
        want_sha,
        have_shas,
        # username=src_username,  # Uncomment if source repo requires auth
        # password=src_password,  # Uncomment if source repo requires auth
    )

    # Upload packfile to destination repository
    send_pack(
        dest_remote_url,
        target_ref,
        from_sha,
        want_sha,
        packfile,
        username=dest_username,
        password=dest_password,
    )


if __name__ == "__main__":
    main()
