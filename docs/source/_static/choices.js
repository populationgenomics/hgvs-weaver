// Split-Pane Reader enhancement for The Choices Weaver Makes (choices.html)
document.addEventListener("DOMContentLoaded", function () {
  const root =
    document.getElementById("the-choices-weaver-makes") ||
    document.getElementById("design-choices");
  if (!root) return;

  document.body.classList.add("choices-split-reader");

  function isDash(val) {
    const t = (val || "").trim();
    return !t || t === "—" || t === "–" || t === "-";
  }

  function statusClass(val) {
    const v = (val || "").toLowerCase();
    if (
      v.indexOf("differs") !== -1 ||
      v.indexOf("rejects") !== -1 ||
      v.indexOf("commits") !== -1 ||
      v.indexOf("errors") !== -1 ||
      v.indexOf("5' in intron") !== -1 ||
      v.indexOf("ext*") !== -1 ||
      v.indexOf("met1?") !== -1 ||
      v.indexOf("spelling") !== -1 ||
      v.indexOf("parses") !== -1 ||
      v.indexOf("text") !== -1
    ) {
      return "v-differ";
    }
    if (v.indexOf("agrees") !== -1) {
      return "v-agree";
    }
    return "v-silent";
  }

  function wrapLeadRuleCallout(sec) {
    const paragraphs = sec.querySelectorAll(":scope > p");
    for (let i = 0; i < paragraphs.length; i++) {
      const p = paragraphs[i];
      const firstChild = p.firstElementChild;
      if (
        firstChild &&
        firstChild.tagName === "STRONG" &&
        p.textContent.trim().indexOf("Example") !== 0
      ) {
        const callout = document.createElement("div");
        callout.className = "choice-rule-callout";
        p.parentNode.insertBefore(callout, p);
        callout.appendChild(p);
        break;
      } else if (p.textContent.trim().indexOf("Example") === 0) {
        break;
      }
    }
    return paragraphs;
  }

  // 1. Parse the "At a glance" table to map #anchor -> {spec, bio, vv, clinvar}
  const glanceMap = {};
  const glanceSection = document.getElementById("at-a-glance");
  if (glanceSection) {
    const rows = glanceSection.querySelectorAll("table tbody tr");
    rows.forEach(function (tr) {
      const cells = tr.querySelectorAll("td");
      if (cells.length >= 5) {
        const link = cells[0].querySelector("a[href^='#']");
        if (link) {
          const anchor = link.getAttribute("href").slice(1);
          glanceMap[anchor] = {
            spec: cells[1].textContent.trim(),
            bio: cells[2].textContent.trim(),
            vv: cells[3].textContent.trim(),
            clinvar: cells[4].textContent.trim(),
          };
        }
        // Colour-code cells in the At-a-glance table as pills
        for (let i = 1; i <= 4; i++) {
          const txt = cells[i].textContent.trim();
          const cls = statusClass(txt);
          if (!isDash(txt)) {
            cells[i].innerHTML =
              '<span class="v-pill ' + cls + '">' + txt + "</span>";
          }
        }
      }
    });
  }

  // 2. Enhance "The principles that decide a choice" cards if present
  const principlesSection = document.getElementById(
    "the-principles-that-decide-a-choice"
  );
  if (principlesSection) {
    const principleCards =
      principlesSection.querySelectorAll(":scope > section");
    const totalPrinciples = principleCards.length;
    principleCards.forEach(function (sec, idx) {
      const h3 = sec.querySelector(":scope > h3");
      if (!h3) return;
      const metaTop = document.createElement("div");
      metaTop.className = "choice-meta-top";
      metaTop.innerHTML =
        '<span class="choice-sec-tag">Principle</span>' +
        '<span class="choice-num-tag">Principle #' +
        (idx + 1) +
        " of " +
        totalPrinciples +
        "</span>";
      sec.insertBefore(metaTop, h3);
      wrapLeadRuleCallout(sec);
    });
  }

  // 3. Collect all h3 choice sections inside h2 domain sections
  const excludedH2 = {
    "at-a-glance": true,
    "how-to-read-this-page": true,
    "the-principles-that-decide-a-choice": true,
  };
  const h2Sections = Array.from(
    root.querySelectorAll(":scope > section")
  ).filter(function (sec) {
    return !excludedH2[sec.id];
  });

  const choiceSections = [];
  h2Sections.forEach(function (h2Sec) {
    const h2 = h2Sec.querySelector(":scope > h2");
    const secName = h2 ? h2.childNodes[0].textContent.trim() : "";
    const subSecs = h2Sec.querySelectorAll(":scope > section");
    subSecs.forEach(function (sub) {
      choiceSections.push({ el: sub, secName: secName });
    });
  });

  const totalChoices = choiceSections.length;

  choiceSections.forEach(function (item, idx) {
    const sec = item.el;
    const h3 = sec.querySelector(":scope > h3");
    if (!h3) return;

    const anchor = sec.id;
    const num = idx + 1;

    // Insert top meta bar before h3
    const metaTop = document.createElement("div");
    metaTop.className = "choice-meta-top";
    metaTop.innerHTML =
      '<span class="choice-sec-tag">' +
      item.secName +
      '</span><span class="choice-num-tag">Choice #' +
      num +
      " of " +
      totalChoices +
      "</span>";
    sec.insertBefore(metaTop, h3);

    // Determine verdicts from At-a-glance or fallback to metadata list
    const g = glanceMap[anchor] || {
      spec: "—",
      bio: "—",
      vv: "—",
      clinvar: "—",
    };

    const allLists = sec.querySelectorAll(":scope > ul");
    const metaList = allLists.length > 0 ? allLists[allLists.length - 1] : null;
    if (metaList) {
      Array.from(metaList.querySelectorAll(":scope > li")).forEach(function (li) {
        const text = li.textContent.trim();
        if (text.indexOf("Spec (") === 0 && isDash(g.spec)) {
          const m = text.match(/^Spec \(([^)]+)\)/);
          if (m) g.spec = m[1];
        } else if (text.indexOf("Tests:") === 0) {
          // Convert Tests bullet into collapsible drawer
          const codeNodes = li.querySelectorAll("code");
          const count = codeNodes.length || 1;
          const details = document.createElement("details");
          details.className = "choice-tests-drawer";
          const innerHtml = li.innerHTML.replace(
            /^<strong>Tests:<\/strong>\s*/,
            ""
          );
          details.innerHTML =
            "<summary>Tests asserting this choice (" +
            count +
            ')</summary><div class="choice-tests-body">' +
            innerHtml +
            "</div>";
          sec.appendChild(details);
          li.remove();
        }
      });
      if (metaList.children.length > 0) {
        metaList.classList.add("choice-meta-list");
      } else {
        metaList.remove();
      }
    }

    // Build verdict pills strip right after h3
    const verdictStrip = document.createElement("div");
    verdictStrip.className = "choice-verdicts";
    const pills = [
      { k: "Spec", v: g.spec },
      { k: "biocommons", v: g.bio },
      { k: "VariantValidator", v: g.vv },
      { k: "ClinVar", v: g.clinvar },
    ];
    pills.forEach(function (p) {
      const cls = statusClass(p.v);
      const span = document.createElement("span");
      span.className = "v-pill " + cls;
      span.innerHTML =
        '<span class="v-k">' + p.k + ":</span> " + (p.v || "—");
      verdictStrip.appendChild(span);
    });
    h3.insertAdjacentElement("afterend", verdictStrip);

    // Wrap the first bold rule paragraph in .choice-rule-callout
    const paragraphs = wrapLeadRuleCallout(sec);

    // Style the "Example..." paragraph and transform text code blocks into Input -> Output tables
    paragraphs.forEach(function (p) {
      if (p.textContent.trim().indexOf("Example") === 0) {
        p.classList.add("choice-example-label");
      }
    });

    const codeBlocks = sec.querySelectorAll(":scope > div.highlight-text");
    codeBlocks.forEach(function (block) {
      const pre = block.querySelector("pre");
      if (!pre) return;
      const lines = pre.textContent.replace(/\n$/, "").split("\n");
      const wrap = document.createElement("div");
      wrap.className = "choice-io-wrap";
      const table = document.createElement("table");
      table.className = "choice-io-table";
      const tbody = document.createElement("tbody");

      lines.forEach(function (line) {
        const tr = document.createElement("tr");
        if (!line.trim()) {
          tr.innerHTML = '<td colspan="4" style="height:6px;border:none;"></td>';
          tbody.appendChild(tr);
          return;
        }
        const arrowUnit = line.indexOf("→") !== -1 ? "→" : null;
        if (arrowUnit) {
          const idx = line.indexOf("→");
          const lhs = line.slice(0, idx).trim();
          const rhsRaw = line.slice(idx + 1).trim();
          const parts = rhsRaw.split(/\s{2,}/);
          const outVal = parts[0] || "";
          const noteVal = parts.slice(1).join("  ");
          const isErr =
            /Error|false|refused|panics|Err\b|not equivalent|Different|rejected/i.test(
              outVal
            );

          const tdIn = document.createElement("td");
          tdIn.className = "choice-io-in";
          tdIn.textContent = lhs;

          const tdArr = document.createElement("td");
          tdArr.className = "choice-io-arr";
          tdArr.textContent = "→";

          const tdOut = document.createElement("td");
          tdOut.className = "choice-io-out" + (isErr ? " is-err" : "");
          tdOut.textContent = outVal;

          const tdNote = document.createElement("td");
          tdNote.className = "choice-io-note";
          tdNote.textContent = noteVal;

          tr.appendChild(tdIn);
          tr.appendChild(tdArr);
          tr.appendChild(tdOut);
          tr.appendChild(tdNote);
        } else {
          const tdRaw = document.createElement("td");
          tdRaw.colSpan = 4;
          tdRaw.className = "choice-io-raw";
          tdRaw.textContent = line;
          tr.appendChild(tdRaw);
        }
        tbody.appendChild(tr);
      });

      table.appendChild(tbody);
      wrap.appendChild(table);
      block.replaceWith(wrap);
    });

    // Decorate matching sidebar navigation item with status dot
    const hasDiff =
      statusClass(g.spec) === "v-differ" ||
      statusClass(g.bio) === "v-differ" ||
      statusClass(g.vv) === "v-differ" ||
      statusClass(g.clinvar) === "v-differ";
    const navLinks = document.querySelectorAll(
      '.wy-menu-vertical a[href$="#' + anchor + '"]'
    );
    navLinks.forEach(function (a) {
      a.classList.add("choice-nav-item");
      const dot = document.createElement("span");
      dot.className =
        "choice-nav-dot " + (hasDiff ? "dot-differ" : "dot-agree");
      dot.title = hasDiff
        ? "Differs from at least one tool or spec"
        : "Agrees across checked tools and spec";
      a.insertBefore(dot, a.firstChild);
    });
  });
});
