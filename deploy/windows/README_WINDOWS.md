# scTCNA web app on Windows Server 2019

The app (`uvicorn app.main:app`) runs from a checkout of this repository as the
Windows service `scTCNA` (NSSM) and listens on `127.0.0.1:8000`. An IIS site forwards
all requests to it and redirects HTTP to HTTPS. The app serves the precomputed tables
in `data\app`; it needs no database. Its memory peaked at 389 MB in the last check
(`docs/memory_check.txt`).

The paths below are examples: `C:\scTCNA` for the repository and `C:\inetpub\sctcna`
for the IIS site folder. Run PowerShell as Administrator.

## You need

| Item | Notes |
| --- | --- |
| Python 3.12, 64-bit, from python.org | install for all users; the pinned numpy has no Windows wheel for 3.11 |
| Git | to clone the repository |
| NSSM pre-release 2.24-101 or newer | [nssm-2.24-101-g897c7ad.zip](https://nssm.cc/ci/nssm-2.24-101-g897c7ad.zip), use `win64\nssm.exe`; the 2.24 release can fail to start services on Windows Server 2019 |
| IIS with URL Rewrite 2.1 and Application Request Routing 3.0 | [URL Rewrite](https://www.iis.net/downloads/microsoft/url-rewrite), [ARR](https://www.iis.net/downloads/microsoft/application-request-routing) |
| A host name and its TLS certificate in the machine store | the site redirects HTTP to HTTPS |

## 1. Get the code

```powershell
git clone https://github.com/BIRL/scTCNA.git C:\scTCNA
```

The tables in `data\app` come with the repository.

## 2. Install

```powershell
cd C:\scTCNA\deploy\windows
powershell -ExecutionPolicy Bypass -File .\setup.ps1
```

`setup.ps1` uses `C:\scTCNA\.venv`, or creates it with Python 3.12, installs
`requirements-app.txt`, and checks that the tables load (`tables load: 10 datasets`).
Options: `-DataDir <folder>` for tables outside the repository, `-Python <path to python.exe>`
when the py launcher does not find Python 3.12.

## 3. Test in the foreground

```powershell
powershell -ExecutionPolicy Bypass -File .\run.ps1
# in a second window:
Invoke-RestMethod http://127.0.0.1:8000/api/health   # status: ok, n_datasets: 10
```

Stop it with Ctrl+C.

## 4. Service (NSSM)

```powershell
$Nssm = 'C:\Tools\nssm-2.24-101-g897c7ad\win64\nssm.exe'
$Root = 'C:\scTCNA'
$Logs = 'C:\ProgramData\scTCNA\logs'
New-Item -ItemType Directory -Force $Logs | Out-Null

& $Nssm install scTCNA "$Root\.venv\Scripts\python.exe" '-m uvicorn app.main:app --host 127.0.0.1 --port 8000 --workers 1'
& $Nssm set scTCNA AppDirectory $Root
& $Nssm set scTCNA AppEnvironmentExtra "SCTCNA_DATA_DIR=$Root\data\app" 'ARROW_DEFAULT_MEMORY_POOL=system' 'PYTHONUNBUFFERED=1' 'PYTHONIOENCODING=utf-8'
& $Nssm set scTCNA AppStdout "$Logs\app.log"
& $Nssm set scTCNA AppStderr "$Logs\app.log"
& $Nssm set scTCNA AppRotateFiles 1
& $Nssm set scTCNA AppRotateBytes 10485760
& $Nssm set scTCNA Start SERVICE_AUTO_START
& $Nssm start scTCNA
```

- `AppDirectory` must be the repository root: the app imports `app` and `sctcna` from there.
- `PYTHONIOENCODING=utf-8`: under the service the output goes to a file instead of a console.
- `AppEnvironmentExtra` replaces the whole list each time it is set.
- `app.log` is rotated when the service starts and the file is over 10 MB.
- The service runs as LocalSystem. To use another account, set it with
  `& $Nssm edit scTCNA` (Log on tab). That account needs read access to `C:\scTCNA` and the
  Python installation, and write access to `C:\ProgramData\scTCNA\logs`.

## 5. IIS site

Once per server, install IIS:

```powershell
Install-WindowsFeature Web-Server -IncludeManagementTools
```

Install URL Rewrite 2.1 and ARR 3.0 (links above), then enable the ARR proxy:

```powershell
& "$env:windir\system32\inetsrv\appcmd.exe" set config -section:system.webServer/proxy /enabled:"True" /commit:apphost
```

Create the site with its own host name, so that it does not clash with Default Web Site on
port 80, and bind the certificate:

```powershell
Import-Module WebAdministration
New-Item -ItemType Directory -Force C:\inetpub\sctcna | Out-Null
Copy-Item C:\scTCNA\deploy\windows\web.config C:\inetpub\sctcna\
New-Website -Name scTCNA -PhysicalPath C:\inetpub\sctcna -Port 80 -HostHeader <host name>
New-WebBinding -Name scTCNA -Protocol https -Port 443 -HostHeader <host name> -SslFlags 1
(Get-WebBinding -Name scTCNA -Protocol https).AddSslCertificate('<certificate thumbprint>', 'My')
```

`web.config` redirects HTTP to HTTPS, forwards every path to `127.0.0.1:8000`, and allows
query strings up to 8192 characters (a 200-gene heatmap URL is longer than the IIS default
of 2048).

## 6. Check

```powershell
Get-Service scTCNA                                   # Running
Invoke-RestMethod http://127.0.0.1:8000/api/health   # status: ok, n_datasets: 10
```

In a browser, open `https://<host name>/api/health`, then the site: Gene, Tumor,
Cross-cancer and Heatmap. On Heatmap, paste 200 genes and press Show, then reload the page
(F5); both must work.

## Updates

```powershell
Stop-Service scTCNA
cd C:\scTCNA
git pull
cd deploy\windows
powershell -ExecutionPolicy Bypass -File .\setup.ps1
Start-Service scTCNA
```

If `web.config` changed, copy it to `C:\inetpub\sctcna` again.

## Troubleshooting

| Symptom | Cause |
| --- | --- |
| `setup.ps1`: "running scripts is disabled" | start it with `powershell -ExecutionPolicy Bypass -File ...` |
| `setup.ps1`: "Stop the service first" | `Stop-Service scTCNA`, then rerun |
| `setup.ps1`: "scTCNA needs 64-bit Python 3.12" for `.venv` | rename `.venv` to `.venv.old` and rerun; `setup.ps1` creates a new environment |
| Service is Paused, or Running while `/api/health` does not answer | read `C:\ProgramData\scTCNA\logs\app.log` and the Application event log (source nssm); usually a wrong `AppDirectory`, or port 8000 already in use |
| `/api/health` answers 503 with "app data not available: ..." | `SCTCNA_DATA_DIR` points to a folder without the app tables (`datasets.parquet` and the other Parquet files), or to tables from an older build; the text after the colon says which |
| IIS error 500.19 | URL Rewrite is not installed |
| IIS error 502.3 | the scTCNA service is not running, or not listening on port 8000 |
| IIS error 404.4 on every path | ARR is not installed or its proxy is not enabled (step 5) |
| Large heatmaps fail with 404.15 through IIS, other pages work | the `web.config` in `C:\inetpub\sctcna` does not set `maxQueryString="8192"`; copy `deploy\windows\web.config` there again |
