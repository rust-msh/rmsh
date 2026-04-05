# Ribbon UI Specification: Implementation Guide

This document provides comprehensive, implementation-level details for building a Ribbon toolbar UI following the Microsoft Office Fluent UI pattern, with specific adaptations for an Ansys AEDT-style engineering application.

---

## 1. Button Types and Visual Specifications

All measurements are at 96 DPI / 100% scaling unless stated otherwise.

### 1.1 Large Button

- **Total cell size**: approximately 44-48 px wide x 82-86 px tall (width varies with label text length)
- **Icon size**: 32x32 px, centered horizontally in the button cell
- **Icon position**: top portion of the button, with ~4 px padding above
- **Label text**: below the icon, centered horizontally
- **Label font**: 8-9pt system font (Segoe UI on Windows)
- **Label wraps**: up to 2 lines maximum; text is center-aligned
- **Padding**: ~4 px horizontal padding on each side of the label; ~2 px between icon and label
- **Minimum width**: ~42 px (for short labels); expands to fit text up to approximately 60 px
- **A large button occupies the full height of the ribbon command area** (all three small-button rows)

### 1.2 Small Button (Medium Variant)

- **Total row height**: ~22 px per row (three rows fit in the ribbon command area)
- **Icon size**: 16x16 px, positioned on the left
- **Label text**: to the right of the icon
- **Padding**: ~3 px between icon and text; ~4 px left/right margins
- **Three small buttons stack vertically** to fill the same height as one large button

### 1.3 Small Button (Compact / Icon-Only Variant)

- **Same row height**: ~22 px
- **Icon size**: 16x16 px
- **No label text** -- command name appears only in the tooltip
- **Used at the "Small" scaling state** when the ribbon is narrowed

### 1.4 Split Button

A composite control combining a primary button (default action) and a secondary dropdown arrow.

- **Primary area**: clicking it executes the default command
- **Dropdown arrow area**: clicking it opens a dropdown menu of related commands
- **Arrow area width**: approximately 12-14 px wide
- **Arrow glyph**: small downward-pointing triangle, approximately 5x3 px, centered in the arrow area
- **Visual separator**: a 1 px line (or subtle color change) divides the primary area from the arrow area
- **Large split button**: the arrow appears below the label text in the bottom portion; the split line is horizontal, separating the icon/label area from the arrow area
- **Small split button**: the arrow appears to the right; the split line is vertical
- **Hover behavior**: the primary area and dropdown area highlight independently on hover, showing which zone will be clicked
- **Use case**: commands with a default plus variants (e.g., Paste / Paste Special, or a Draw Primitive with options)

### 1.5 Toggle Button

Functionally identical to a regular button but maintains a checked/unchecked state.

- **Visual size**: same as the equivalent button size (large or small)
- **Checked/active state**: persistent highlighted background with a visible border (see Button States section)
- **Checked + hover**: a slightly different shade to distinguish from plain checked
- **Use case**: Bold, Italic, View toggles, grid display on/off

### 1.6 Dropdown Button

The entire button opens a dropdown menu -- there is no separate primary action.

- **Visually similar to a split button** but without the visual separator between button face and arrow
- **The dropdown arrow glyph is part of the button face** itself, typically to the right of the label (small) or below the label (large)
- **Every click on any part of the button opens the menu**
- **Use case**: command lists with no single default (e.g., "New from template", "Insert Component")

### 1.7 Gallery (In-Ribbon Gallery)

A grid of visual thumbnail options displayed directly within the ribbon or as a dropdown.

- **In-Ribbon Gallery**: a horizontal strip of thumbnail items shown inline in the ribbon command area, with up/down scroll arrows on the right side and a "More" dropdown arrow below them
- **Gallery item sizes**: varies by content -- common sizes include 48x48, 72x96, 96x72, 96x96 px
- **Dropdown gallery**: when the "More" arrow is clicked, a larger popup shows the full grid of options
- **Gallery popup**: may include category headers (non-clickable text labels), separator lines, and additional commands at the bottom
- **Gallery scaling**: in narrow ribbon states, the in-ribbon gallery collapses to a single dropdown button showing just the currently selected item's thumbnail
- **AEDT context**: material selection, color/fill selection, plot style selection

### 1.8 Combo Box / Text Input in Ribbon

- **Combo Box**: dropdown list with an editable text field, using the "input family" control type
- **Height**: matches small-button row height (~22 px)
- **Width**: typically 100-150 px for the text area plus ~16 px for the dropdown arrow
- **Border**: 1 px solid border in normal state; highlighted border on focus
- **Use case in AEDT**: unit selection, coordinate entry, design parameter values
- **Spinner**: numeric input with up/down increment buttons, also part of the input family

### 1.9 Checkbox in Ribbon

- **Yes, the Windows Ribbon Framework supports CheckBox controls directly in the ribbon**
- **The `ThreeButtonsAndOneCheckBox` SizeDefinition template** explicitly includes a CheckBox alongside three buttons
- **Appearance**: standard checkbox with a text label, fitting within the small-button row height
- **Use case**: enabling/disabling features inline (e.g., "Show Grid", "Snap to Grid")

---

## 2. Button States

### 2.1 Normal (Rest / Default)

- **Background**: transparent or matches the ribbon background (typically a light gray, e.g., #F5F5F5)
- **Border**: none / invisible
- **Icon**: full color, full opacity
- **Text**: standard color (dark gray/black, e.g., #333333)
- **Appearance**: flat, minimal chrome -- the button blends into the ribbon surface

### 2.2 Hover (Mouse Over)

- **Background**: subtle highlight fill appears -- typically a light warm color (Office uses a pale blue or pale orange tint, e.g., #E5F1FB or #FDE8D0)
- **Border**: 1 px solid border appears around the button (slightly darker than the fill, e.g., #B8D6FB)
- **Icon**: unchanged (full color)
- **Text**: unchanged
- **Transition**: instantaneous or very fast (<50ms)
- **Split button**: only the hovered zone (primary or arrow) highlights; the other zone remains in rest state

### 2.3 Pressed (Mouse Down / Active)

- **Background**: darker fill than hover (e.g., #CCE4F7 or a medium blue tint)
- **Border**: 1 px solid, slightly darker than hover border
- **Slight visual "inset" effect** -- the background color is noticeably darker, giving a pressed/depressed feel
- **Icon and text**: unchanged
- **Duration**: only while mouse button is held down

### 2.4 Disabled / Grayed

- **Icon**: grayscale or reduced opacity (~40-50% opacity)
- **Text**: lighter color (e.g., #AAAAAA or ~50% of normal text opacity)
- **Background**: no highlight on hover (hover effects are suppressed)
- **Cursor**: default arrow (not the hand/pointer; some implementations use `not-allowed`)
- **Border**: never appears
- **No interaction**: clicks do not fire, no tooltip (or a simpler "not available" tooltip)

### 2.5 Checked / Active (Toggle Buttons)

- **Background**: persistent highlight fill -- typically a saturated accent color (light blue, e.g., #D0E8FF, or in Office 2007-style a light orange)
- **Border**: 1 px solid accent-colored border, always visible
- **Icon and text**: full color, unchanged
- **Checked + Hover**: slightly different shade than checked-rest (e.g., a bit darker or warmer) to indicate interactive responsiveness
- **Checked + Pressed**: even darker shade momentarily

### 2.6 Focused (Keyboard Focus)

- **Dotted or dashed focus rectangle** around the button, or a visible focus ring
- **Background**: may show a subtle highlight similar to hover
- **Keytip overlay**: when Alt is pressed, keytip badges appear over buttons showing keyboard shortcut letters
- **Tab navigation**: moves focus between controls; Enter/Space activates the focused control

---

## 3. Group Layout Rules

### 3.1 Ribbon Anatomy (Overall Dimensions at 96 DPI)

```
+----------------------------------------------------------+
| Tab Strip (25-30 px tall)                                |
| [Desktop] [Draw] [Model] [Simulation] [Results] [View]  |
+----------------------------------------------------------+
| Ribbon Command Area (~72-80 px tall)                     |
|                                                          |
| +--------+ +---------+ | +----+ +----+ +----+ |         |
| | Large  | | Large   | | |Sm 1| |Sm 4| |Sm 7| |         |
| | Button | | Button  | | |Sm 2| |Sm 5| |Sm 8| |         |
| |  32x32 | |  32x32  | | |Sm 3| |Sm 6| |Sm 9| |         |
| | Label  | | Label   | | +----+ +----+ +----+ |         |
| +--------+ +---------+ | Group Separator (1px)|         |
|                         |                      |         |
| [-- Group Label --]     [--- Group Label ---]  |         |
+----------------------------------------------------------+
| Group Label Strip (~15-18 px tall)                       |
+----------------------------------------------------------+
```

**Total ribbon height** (tabs + commands + group labels): approximately **115-130 px** in modern Office; the original 2007 spec was ~91-100 px.

### 3.2 Standard Layout Patterns

**Pattern A: Single Large Button (OneButton template)**
- One large button-family control occupying the full group width
- Only Large size is supported (no scaling)

**Pattern B: Two or Three Large Buttons Side by Side (TwoButtons, ThreeButtons)**
- 2-3 large buttons arranged horizontally
- At Medium scaling: buttons shrink to small icons with text labels in stacked rows

**Pattern C: One Large + Stack of Small (ThreeButtons-OneBigAndTwoSmall)**
- The first button stays large (prominent)
- The remaining two display as small stacked buttons to its right
- Supports Large, Medium, and Small group sizes

**Pattern D: Stack of 3 Small Buttons (Three-Row Layout)**
- Three small buttons stacked vertically (one per row)
- Each row is ~22 px tall
- This is the standard layout when groups are at Medium size

**Pattern E: Mixed Arrangements (FourButtons through ElevenButtons)**
- Microsoft provides predefined templates for 4-11 buttons
- At Large size: mix of large and small buttons
- At Medium size: stacked rows of small buttons with labels
- At Small size: stacked rows of small icon-only buttons

**Pattern F: ButtonGroups (Complex Grid)**
- Up to 32 button-family controls arranged in two rows of control groups
- Supports an optional full-size button alongside dense grids of small buttons
- Used for formatting toolbar-style layouts (like Office's paragraph/font sections)

**Pattern G: BigButtonsAndSmallButtonsOrInputs**
- Up to 2 large buttons followed by 2-3 small buttons or input controls (ComboBox, Spinner)
- Only Large and Medium sizes supported

**Pattern H: Gallery + Buttons (InRibbonGalleryAndButtons-GalleryScalesFirst)**
- One in-ribbon gallery plus 2-3 button-family controls
- The gallery collapses to a popup first, before the buttons shrink

### 3.3 Group Separator

- **Style**: 1 px vertical line
- **Color**: a medium gray, slightly darker than the ribbon background (e.g., #C0C0C0 or themed)
- **Height**: spans the full height of the ribbon command area (not into the group label strip)
- **Padding**: approximately 3-4 px on each side of the separator line
- **Total separator space**: ~6-8 px wide

### 3.4 Group Label

- **Position**: bottom of the group, within the group label strip
- **Font size**: ~8pt or 11px, Segoe UI or system font
- **Alignment**: center-aligned horizontally within the group
- **Background**: slightly different shade than the command area (subtle distinction), or transparent
- **Text color**: medium gray (#666666)
- **Height**: ~15-18 px
- **An optional dropdown arrow** may appear in the group label area (called a "dialog launcher") to open a related dialog

---

## 4. Dropdown Menu Details

### 4.1 Menu Appearance

- **Border**: 1 px solid (#D1D1D1 or themed)
- **Border radius**: 4 px (Fluent 2 standard; older Office used 0-2 px)
- **Background**: white (#FFFFFF) or theme surface color
- **Shadow**: `box-shadow: 0px 4px 8px rgba(0,0,0,0.14), 0px 0px 2px rgba(0,0,0,0.12)` (Fluent 2 `$shadow8`)
- **Max width**: ~300 px for standard command menus; galleries can be wider
- **Appears**: directly below the button that triggered it, aligned to the left edge

### 4.2 Menu Item

- **Height**: 32 px (standard) or 24 px (compact/dense mode)
- **Padding**: 0 12 px horizontal; vertically centered content
- **Structure**: [16x16 icon] [12px gap] [label text] [flexible space] [shortcut key text] [optional submenu arrow]
- **Icon area**: 16x16 px, with ~4 px margin on the left
- **Shortcut key text**: right-aligned, lighter color (#888888), font size same as label
- **Submenu arrow**: small right-pointing triangle (~5x8 px) on the far right
- **Text**: left-aligned, standard font (Segoe UI 9pt)
- **Hover**: background highlight (same light blue/warm tint as button hover), covers the full width of the item

### 4.3 Menu Separator

- **Horizontal line**: 1 px solid (#E0E0E0)
- **Margin**: ~4 px above and below the line
- **Indent**: may be indented to align with the text (past the icon area), or full width

### 4.4 Menu Header / Category Label

- **Non-clickable text**, typically bold or semi-bold
- **Background**: may have a slightly tinted background
- **Height**: same as a menu item (~32 px)
- **No hover effect**
- **Used to label groups of related commands within a single menu**

### 4.5 Submenu Behavior

- **Hover delay**: approximately 200-400 ms before the submenu appears (prevents accidental triggering)
- **Position**: the submenu appears to the right of the parent menu, aligned to the top of the hovered item
- **If insufficient space on the right**: the submenu appears to the left instead
- **Moving the mouse diagonally**: a brief "tolerance zone" allows moving toward the submenu without the selection jumping to adjacent items

### 4.6 Gallery Dropdown

- **Grid layout**: items arranged in a grid (e.g., 4-8 columns)
- **Item size**: larger thumbnails (48x48, 64x64, or custom sizes)
- **Category headers**: text labels above groups of gallery items
- **Scroll**: if the gallery is large, a vertical scrollbar appears
- **Bottom command area**: additional buttons below the grid (e.g., "More Options...", "Browse...")

### 4.7 Recent Items List (Desktop/File Tab)

- **Item height**: taller than standard menu items (~40-48 px)
- **Structure**: [file icon] [file name (bold)] [file path (lighter, smaller text)]
- **May show two lines**: name on first line, path on second
- **Pin icon**: on the right side to pin/unpin frequently used items
- **Typically shown in the Backstage/File area**, not in a standard dropdown

---

## 5. Ribbon Tab Behavior

### 5.1 Tab Strip

- **Height**: approximately 25-30 px
- **Font**: Segoe UI 9pt or equivalent system font
- **Text color (inactive tab)**: medium gray (#444444)
- **Text color (active tab)**: dark / black
- **Active tab indicator**: the active tab has a bottom border that merges with the ribbon body below (same background color), creating a connected appearance. In some implementations, the active tab has a colored accent bar on top (2-3 px).

### 5.2 Tab Hover

- **Background**: subtle highlight on hover (lighter than the active tab)
- **Transition**: smooth, ~100ms
- **Underline/border**: a subtle bottom-border highlight may appear

### 5.3 Active Tab Connection

- **The active tab's bottom edge has no bottom border** (or shares the same background as the ribbon command area), making it visually "connected" to the content panel below
- **Inactive tabs** have a visible bottom border separating them from the ribbon body
- **This is a key visual detail** -- the selected tab appears to be physically part of the ribbon panel

### 5.4 Contextual Tabs

- **Contextual tabs appear only when relevant** (e.g., selecting an object shows a "Format" tab)
- **Colored accent bar**: a 2-3 px colored bar above the contextual tab group (orange, green, blue, etc.) to distinguish them from permanent tabs
- **Tab text may be colored** to match the accent
- **A group label above the contextual tabs** shows the context category name (e.g., "Drawing Tools", "Table Tools")
- **AEDT behavior**: the Draw, Model, Results tabs are design-type-contextual -- they change based on the active solver (HFSS vs Maxwell vs Circuit). In AEDT this is more like persistent context switching than transient contextual tabs.

---

## 6. Adaptive Collapse Behavior (Responsive Sizing)

The ribbon uses a defined **ScalingPolicy** that specifies the exact order in which groups reduce in size.

### 6.1 Four Size States

| State | Icon | Label | Description |
|-------|------|-------|-------------|
| **Large** | 32x32 | Below icon | Full-size presentation; primary commands |
| **Medium** | 16x16 | Beside icon | Compact rows; stacked up to 3 per column |
| **Small** | 16x16 | Hidden | Icon-only; command name in tooltip only |
| **Popup** | N/A | N/A | Entire group collapses to a single dropdown button |

### 6.2 Collapse Sequence

As the window narrows, groups shrink in a developer-specified priority order:

1. **Stage 1 (Large -> Medium)**: lower-priority groups shrink first. Large buttons become small-with-label rows. Gallery controls may lose columns.
2. **Stage 2 (Medium -> Small)**: labels are hidden, leaving icon-only buttons in stacked rows.
3. **Stage 3 (Small -> Popup)**: the entire group collapses into a single popup button.

The order is specified per-group. Typically:
- **Most-used groups stay Large longest** (e.g., Clipboard, Primary Draw tools)
- **Less-used groups collapse first** (e.g., View options, Formatting details)
- **Each step down is explicitly declared** in the ScalingPolicy, not automatic

### 6.3 Collapsed Group (Popup) Appearance

- **Single button** approximately 48-56 px wide x full ribbon command area height
- **Displays**: the group's icon (32x32) + the group name text below
- **A small dropdown arrow** below the label
- **On click**: opens a popup panel showing the group's controls in their Large layout
- **The popup panel has the same appearance as the ribbon command area** (same background, same control styling)

### 6.4 ScalingPolicy Example (XML)

```xml
<ScalingPolicy>
  <ScalingPolicy.IdealSizes>
    <Scale Group="GroupClipboard" Size="Large"/>
    <Scale Group="GroupFont" Size="Large"/>
    <Scale Group="GroupParagraph" Size="Large"/>
  </ScalingPolicy.IdealSizes>
  <!-- Paragraph shrinks first -->
  <Scale Group="GroupParagraph" Size="Medium"/>
  <Scale Group="GroupFont" Size="Medium"/>
  <Scale Group="GroupParagraph" Size="Small"/>
  <Scale Group="GroupFont" Size="Small"/>
  <Scale Group="GroupParagraph" Size="Popup"/>
  <Scale Group="GroupFont" Size="Popup"/>
</ScalingPolicy>
```

---

## 7. Quick Access Toolbar (QAT)

### 7.1 Position

- **Default**: above the ribbon, in the title bar area
- **Optional**: can be placed below the ribbon (user preference)
- **Height**: matches the title bar height (~20-24 px)

### 7.2 Button Style

- **Icon-only**: 16x16 px icons
- **Hit target area**: approximately 22x22 px per button (including padding)
- **No labels**: command names appear in tooltips
- **Background**: transparent (blends with title bar)
- **Spacing between buttons**: ~2-3 px

### 7.3 Customization Dropdown

- **Located at the right end of the QAT**
- **Small dropdown arrow** (~8 px wide area)
- **Opens a menu** with: commonly added commands, "More Commands...", "Show Below the Ribbon", etc.

### 7.4 Separator

- **Thin vertical line** (1 px, light gray)
- **Padding**: ~2-4 px on each side
- **Total space**: ~6-9 px
- **Used to visually group related QAT commands**

---

## 8. Tooltip / Screentip Behavior

### 8.1 Simple Tooltip

- **Content**: just the command name (e.g., "Paste")
- **Appears after**: ~500 ms hover delay (system default initial delay)
- **Position**: below the button
- **Duration**: stays visible for ~5000 ms (autopop duration), then fades
- **Reshow delay**: when moving between controls with tooltip already visible, next tooltip appears after ~100 ms

### 8.2 Enhanced Screentip (Rich Tooltip)

Office and AEDT-style applications support richer tooltips with:

- **Title**: bold command name at the top (e.g., "Paste (Ctrl+V)")
- **Description**: 1-3 lines of explanatory text below the title
- **Optional image**: a preview or illustration (~100-200 px wide)
- **Optional shortcut key**: shown after the title or on a separate line
- **Separator line**: between the title and description
- **"Tell me more" / help link**: at the bottom, linking to documentation

### 8.3 Tooltip Timing (Win32 Platform Standard)

| Phase | Default Value | Description |
|-------|--------------|-------------|
| **Initial delay** | 500 ms | Time cursor must hover before tooltip appears |
| **Autopop (show duration)** | 5000 ms | How long the tooltip stays visible |
| **Reshow delay** | 100 ms | Delay when moving between controls while tooltips are active |

These values are based on the system double-click time and scale proportionally.

### 8.4 AEDT Tooltips

- AEDT uses enhanced screentips similar to Office
- **Command name** is shown in bold as the title
- **Brief description** explains what the command does
- **Shortcut key** (if any) is displayed
- The tooltip appears below the ribbon button, positioned to avoid going off-screen

---

## 9. AEDT-Specific Ribbon Details

### 9.1 Tab Structure

AEDT uses a combination of permanent and context-dependent tabs:

**Permanent Tabs (always visible):**

| Tab | Purpose |
|-----|---------|
| Desktop | Project operations: new/open/save, insert design (HFSS, Maxwell, Q3D, Circuit, etc.), import/export |
| View | Zoom, pan, rotate, visibility settings, window layout, coordinate system display (Large/Small/Hide triad) |
| Simulation | Analysis setup, frequency sweeps, HPC configuration, Validate, Analyze All |
| Automation | ACT toolkits, scripting, macro record/playback |

**Context-Dependent Tabs (change based on active solver):**

| Solver | Contextual Tabs |
|--------|----------------|
| HFSS, Maxwell, Q3D, Icepak, Mechanical | Draw, Model, Results |
| HFSS 3D Layout | Layout, Results |
| Circuit / Twin Builder | Schematic, Results |

### 9.2 Draw Tab (HFSS/Maxwell/Q3D)

Groups and controls:

**Primitives Group (Large buttons):**
- Box, Cylinder, Sphere, Cone, Torus, Helix, Regular Polyhedron
- These are the most-used commands and appear as large buttons with 32x32 icons

**Lines/Curves Group (Small/Medium buttons):**
- Line, Arc (3-point, center-point), Spline, Circle, Ellipse, Rectangle, Regular Polygon
- Typically small buttons stacked in rows

**Boolean Operations Group (Small buttons):**
- Unite, Subtract, Intersect, Split
- Small stacked buttons

**Edit Geometry Group:**
- Sweep Along Vector, Sweep Along Path, Sweep Around Axis
- Section, Clone, Mirror, Move, Rotate, Scale, Offset

**Coordinate System Group:**
- Create Relative CS, Create Face CS
- May include a dropdown for coordinate system type selection

### 9.3 Model Tab (HFSS/Maxwell/Q3D)

**Boundaries Group:**
- Assign Boundary (dropdown with submenu of boundary types)
- For HFSS: Perfect E, Perfect H, Impedance, Radiation, PML, etc.
- For Maxwell: various boundary condition types

**Excitations Group:**
- Assign Excitation (dropdown with submenu)
- For HFSS: Wave Port, Lumped Port, Floquet Port, etc.

**Mesh Operations Group:**
- Assign Mesh Operation (dropdown)
- Length-Based, Skin Depth, Curvilinear, etc.

**Materials Group:**
- Assign Material (may include a gallery or dropdown with material list)

**Design Settings:**
- Solution Type, Design Properties, Parametric Setup

### 9.4 Simulation Tab

**Setup Group:**
- Add Solution Setup (large button or split button)
- Add Frequency Sweep
- Edit Setup / Edit Sweep

**Validation Group:**
- Validate (large button with checkmark icon)

**Analysis Group:**
- Analyze All (large button, prominent -- primary action)

**HPC Group:**
- HPC and Analysis Options
- Job Monitor

### 9.5 Results Tab

**Reports Group (Large buttons):**
- Create Report (dropdown with report types: Rectangular Plot, Polar Plot, Smith Chart, Data Table, 3D Polar Plot)
- Solution Data

**Field Visualization Group:**
- Field Overlays (dropdown: plot E-field, H-field, Current Density, SAR, etc.)
- Far-Field Plot (2D, 3D radiation patterns)

**Animation Group:**
- Animate (create field animations)

**Export Group:**
- Export results to file

### 9.6 Desktop Tab

**Design Insertion Group (Large buttons with dropdowns):**
- HFSS (insert new HFSS design)
- Maxwell 2D, Maxwell 3D
- Q3D Extractor
- Circuit Design
- Twin Builder (Simplorer)
- Icepak
- Mechanical
- Each of these may be a dropdown button offering additional options (e.g., "Insert from Library")

**Project Group:**
- New, Open, Save, Save As
- Import, Export

### 9.7 Special Controls in AEDT

**Status Bar Controls (bottom of the modeler window, not in the ribbon):**
- X, Y, Z coordinate entry text boxes
- Absolute / Relative toggle selector
- Cartesian / Cylindrical / Spherical coordinate system selector
- Unit selector dropdown (mm, cm, mil, in, etc.)

**These are NOT ribbon controls** but are important UI elements in AEDT's geometry entry workflow.

**Unit selector in ribbon:**
- A ComboBox or dropdown in the Draw tab's group
- Allows changing the model's working units

**Material selector:**
- May use a gallery-style dropdown showing material names with color swatches
- Or a standard dropdown opening a material browser dialog

---

## 10. Color Specifications Reference

### 10.1 Ribbon Chrome Colors (Light Theme, Office-Style)

| Element | Color |
|---------|-------|
| Ribbon background | #F5F5F5 (very light gray) |
| Tab strip background | #FFFFFF or slightly lighter than ribbon |
| Active tab background | same as ribbon background (connected) |
| Inactive tab text | #444444 |
| Active tab text | #000000 or #1A1A1A |
| Group label text | #666666 |
| Group separator line | #C0C0C0 |
| Button hover fill | #E5F1FB (pale blue) |
| Button hover border | #B8D6FB |
| Button pressed fill | #CCE4F7 |
| Button pressed border | #98C4EA |
| Toggle checked fill | #D0E8FF |
| Toggle checked border | #6CA0DC |
| Disabled text | #AAAAAA |
| Disabled icon opacity | 40-50% |
| Menu background | #FFFFFF |
| Menu border | #D1D1D1 |
| Menu item hover | #E5F1FB |
| Menu separator | #E0E0E0 |

### 10.2 Dark Theme Adaptations (if needed)

For a dark theme, invert the luminance values:
- Ribbon background: #2D2D2D
- Text: #E0E0E0
- Hover fill: #3D3D3D with accent tint
- Borders: #555555
- Active tab: #2D2D2D (matches ribbon)
- Menu background: #333333

---

## 11. Implementation Notes for CSS/Web Framework

### 11.1 Key CSS Variables

```css
:root {
  /* Ribbon chrome */
  --ribbon-bg: #F5F5F5;
  --ribbon-tab-height: 28px;
  --ribbon-command-height: 78px;
  --ribbon-group-label-height: 16px;
  --ribbon-total-height: calc(var(--ribbon-tab-height) + var(--ribbon-command-height) + var(--ribbon-group-label-height));

  /* Icons */
  --icon-large: 32px;
  --icon-small: 16px;

  /* Button sizing */
  --btn-large-width: 48px;
  --btn-large-height: var(--ribbon-command-height);
  --btn-small-row-height: 22px;
  --btn-padding-h: 4px;

  /* Group separator */
  --group-separator-width: 1px;
  --group-separator-margin: 4px;

  /* Interactive states */
  --hover-bg: #E5F1FB;
  --hover-border: #B8D6FB;
  --pressed-bg: #CCE4F7;
  --pressed-border: #98C4EA;
  --checked-bg: #D0E8FF;
  --checked-border: #6CA0DC;
  --disabled-opacity: 0.45;

  /* Menu */
  --menu-item-height: 32px;
  --menu-border-radius: 4px;
  --menu-shadow: 0px 4px 8px rgba(0,0,0,0.14), 0px 0px 2px rgba(0,0,0,0.12);

  /* Split button */
  --split-arrow-width: 14px;
  --split-separator: 1px solid #D0D0D0;

  /* Tooltip */
  --tooltip-delay: 500ms;
  --tooltip-duration: 5000ms;
}
```

### 11.2 Predefined SizeDefinition Templates (from Windows Ribbon Framework)

These templates define how controls are arranged within a group at each size state:

| Template Name | Controls | Supported Sizes |
|--------------|----------|-----------------|
| OneButton | 1 button | Large only |
| TwoButtons | 2 buttons | Large, Medium |
| ThreeButtons | 3 buttons | Large, Medium |
| ThreeButtons-OneBigAndTwoSmall | 3 buttons (1st prominent) | Large, Medium, Small |
| ThreeButtonsAndOneCheckBox | 3 buttons + 1 checkbox | Large, Medium |
| FourButtons | 4 buttons | Large, Medium, Small |
| FiveButtons | 5 buttons | Large, Medium, Small |
| FiveOrSixButtons | 5-6 buttons | Large, Medium, Small |
| SixButtons | 6 buttons | Large, Medium, Small |
| SixButtons-TwoColumns | 6 buttons (alt layout) | Large, Medium, Small |
| SevenButtons | 7 buttons | Large, Medium, Small |
| EightButtons | 8 buttons | Large, Medium, Small |
| EightButtons-LastThreeSmall | 8 buttons (last 3 grouped) | Large, Medium, Small |
| NineButtons | 9 buttons | Large, Medium, Small |
| TenButtons | 10 buttons | Large, Medium, Small |
| ElevenButtons | 11 buttons | Large, Medium, Small |
| OneFontControl | 1 FontControl | Large, Medium |
| OneInRibbonGallery | 1 InRibbonGallery | Large, Small |
| InRibbonGalleryAndBigButton | 1 gallery + 1 button | Large, Small |
| InRibbonGalleryAndButtons-GalleryScalesFirst | 1 gallery + 2-3 buttons | Large, Medium, Small |
| ButtonGroups | Up to 32 buttons in grid | Large, Medium, Small |
| ButtonGroupsAndInputs | 2 inputs + 29 buttons | Large, Medium |
| BigButtonsAndSmallButtonsOrInputs | 2 big + 2-3 small/input | Large, Medium |

**Button family controls**: Button, Toggle Button, Drop-Down Button, Split Button, Drop-Down Gallery, Split Button Gallery, Drop-Down Color Picker.

**Input family controls**: Combo Box, Spinner.

**Standalone controls**: CheckBox, In-Ribbon Gallery (used only where explicitly allowed in a template).

---

## 12. Sources and References

- [Windows Ribbon Framework: Size Definitions and Scaling Policies](https://learn.microsoft.com/en-us/windows/win32/windowsribbon/windowsribbon-templates)
- [Windows Ribbon Framework: ScalingPolicy Element](https://learn.microsoft.com/en-us/windows/win32/windowsribbon/windowsribbon-element-scalingpolicy)
- [Windows Ribbon Framework: Scale Element](https://learn.microsoft.com/en-us/windows/win32/windowsribbon/windowsribbon-element-scale)
- [Windows Ribbon Framework: Split Button Control](https://learn.microsoft.com/en-us/windows/win32/windowsribbon/windowsribbon-controls-splitbutton)
- [Windows Ribbon Framework: Toggle Button Control](https://learn.microsoft.com/en-us/windows/win32/windowsribbon/windowsribbon-controls-togglebutton)
- [Windows Ribbon Framework: Button Element](https://learn.microsoft.com/en-us/windows/win32/windowsribbon/windowsribbon-element-button)
- [Windows Ribbon Framework: GroupSizeDefinition](https://learn.microsoft.com/en-us/windows/win32/windowsribbon/windowsribbon-element-groupsizedefinition)
- [Office Add-in Icon Guidelines](https://learn.microsoft.com/en-us/office/dev/add-ins/design/add-in-icons-fresh)
- [Fluent 2 Design System: Elevation (Shadows)](https://fluent2.microsoft.design/elevation)
- [Fluent 2 Design System: Shapes (Border Radius)](https://fluent2.microsoft.design/shapes)
- [Fluent 2 Design System: React Menu Component](https://fluent2.microsoft.design/components/web/react/core/menu/usage)
- [Win32 Tooltip Timing: TTM_SETDELAYTIME](https://learn.microsoft.com/en-us/windows/win32/controls/ttm-setdelaytime)
- [Win32 UX Guide: Tooltips and Infotips](https://learn.microsoft.com/en-us/windows/win32/uxguide/ctrl-tooltips-and-infotips)
- [Actipro WPF Ribbon UI Guidelines](https://www.actiprosoftware.com/docs/controls/wpf/ribbon/ribbonui-guidelines)
- [Fluent.Ribbon Sizing Concepts](https://fluentribbon.github.io/documentation/concepts/sizing)
- [Metro UI 5.1 Ribbon Menu (CSS Implementation)](https://v5.metroui.org.ua/components/ribbon-menu/)
- [Ansys AEDT: Working with Ribbons (HFSS)](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v242/en/Subsystems/HFSS/Content/GettingStarted/WorkingWithRibbons.htm)
- [Ansys AEDT: Working with Ribbons (Icepak)](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v242/en/Subsystems/Icepak/Content/GettingStarted/WorkingWithRibbons.htm)
- [Ansys Innovation Courses: Intro to AEDT User Interface](https://innovationspace.ansys.com/courses/courses/intro-to-ansys-hfss/lessons/intro-to-aedt-user-interface-lesson-1/)
- [2007 Microsoft Office Fluent UI Design Guidelines (archived)](https://www.scribd.com/document/53387618/2007-Microsoft-Office-Fluent-UI-Design-Guidelines-License)
- [Icon Dimensions in the Office Ribbon (nolongerset.com)](https://nolongerset.com/icon-dimensions-in-the-office-ribbon/)
