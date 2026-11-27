; NiTCAD Windows installer (Inno Setup 6). Built by packaging/build-installer.cmd, which stages
; `cmake --install` into build\stage and passes the version: iscc /DAppVersion=0.1.0 nitcad.iss
#ifndef AppVersion
  #define AppVersion "0.1.0"
#endif

[Setup]
AppId={{5C7E2B1A-8F3D-4E6A-9B2C-1D4F7A9E3C60}
AppName=NiTCAD
AppVersion={#AppVersion}
AppPublisher=NiTCAD
DefaultDirName={autopf}\NiTCAD
DefaultGroupName=NiTCAD
; Per-user by default (no admin prompt); the user may choose all users.
PrivilegesRequired=lowest
PrivilegesRequiredOverridesAllowed=dialog
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
; Direct3D 12 and the UCRT: Windows 10 or later.
MinVersion=10.0
LicenseFile=..\LICENSE
SetupIconFile=..\app\icon\nitcad.ico
UninstallDisplayIcon={app}\nitcad-app.exe
OutputDir=..\build\installer
OutputBaseFilename=NiTCAD-{#AppVersion}-setup
Compression=lzma2
SolidCompression=yes
WizardStyle=modern
ChangesEnvironment=yes

[Tasks]
Name: "desktopicon"; Description: "Create a &desktop shortcut"; GroupDescription: "Shortcuts:"
Name: "addtopath"; Description: "Add the nitcad command-line tool to PATH"; GroupDescription: "Command line:"; Flags: unchecked

[Files]
Source: "..\build\stage\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
Name: "{group}\NiTCAD"; Filename: "{app}\nitcad-app.exe"
Name: "{group}\Uninstall NiTCAD"; Filename: "{uninstallexe}"
Name: "{autodesktop}\NiTCAD"; Filename: "{app}\nitcad-app.exe"; Tasks: desktopicon

[Registry]
Root: HKA; Subkey: "{code:EnvKey}"; ValueType: expandsz; ValueName: "Path"; ValueData: "{olddata};{app}"; Tasks: addtopath; Check: NeedsPath(ExpandConstant('{app}'))

[Run]
Filename: "{app}\nitcad-app.exe"; Description: "Launch NiTCAD"; Flags: nowait postinstall skipifsilent

[Code]
function EnvKey(Param: String): String;
begin
  if IsAdminInstallMode then
    Result := 'SYSTEM\CurrentControlSet\Control\Session Manager\Environment'
  else
    Result := 'Environment';
end;

function RootKey: Integer;
begin
  if IsAdminInstallMode then Result := HKEY_LOCAL_MACHINE else Result := HKEY_CURRENT_USER;
end;

function NeedsPath(Dir: String): Boolean;
var
  Old: String;
begin
  if not RegQueryStringValue(RootKey, EnvKey(''), 'Path', Old) then
  begin
    Result := True;
    exit;
  end;
  Result := Pos(';' + Uppercase(Dir) + ';', ';' + Uppercase(Old) + ';') = 0;
end;

procedure CurUninstallStepChanged(CurUninstallStep: TUninstallStep);
var
  Old, Dir: String;
  P: Integer;
begin
  if CurUninstallStep <> usPostUninstall then exit;
  if not RegQueryStringValue(RootKey, EnvKey(''), 'Path', Old) then exit;
  Dir := ExpandConstant('{app}');
  P := Pos(';' + Uppercase(Dir), Uppercase(Old));
  if P > 0 then
  begin
    Delete(Old, P, Length(Dir) + 1);
    RegWriteExpandStringValue(RootKey, EnvKey(''), 'Path', Old);
  end;
end;
